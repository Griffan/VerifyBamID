#include "SimplePileupViewer.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>
#include <errno.h>
#include <algorithm>
#include <sys/stat.h>
#include <getopt.h>
#include "htslib/sam.h"

#include "htslib/kstring.h"
#include "htslib/khash_str2int.h"
#include "sam_header.h"
#include "samtools.h"
#include <assert.h>
#include <sample.h>
#include <Error.h>
#include "htslib/vcf.h"
#include "sample.h"


static inline int printw(int c, FILE *fp) {
    char buf[16];
    int l, x;
    if (c == 0) return fputc('0', fp);
    for (l = 0, x = c < 0 ? -c : c; x > 0; x /= 10) buf[l++] = x % 10 + '0';
    if (c < 0) buf[l++] = '-';
    buf[l] = 0;
    for (x = 0; x < l / 2; ++x) {
        int y = buf[x];
        buf[x] = buf[l - 1 - x];
        buf[l - 1 - x] = y;
    }
    fputs(buf, fp);
    return 0;
}

static inline void
pileup_seq(FILE *fp, const bam_pileup1_t *p, int pos, int ref_len, const char *ref, std::vector<char> &tmpBase) {
    int j;
    if (p->is_head) {
//        putc('^', fp);
//        putc(p->b->core.qual > 93? 126 : p->b->core.qual + 33, fp);
    }
    if (!p->is_del) {
        int c = p->qpos < p->b->core.l_qseq
                ? seq_nt16_str[bam_seqi(bam_get_seq(p->b), p->qpos)]
                : 'N';
        if (ref) {
            int rb = pos < ref_len ? ref[pos] : 'N';
            if (c == '=' || seq_nt16_table[c] == seq_nt16_table[rb]) c = bam_is_rev(p->b) ? ',' : '.';
            else c = bam_is_rev(p->b) ? tolower(c) : toupper(c);
        } else {
            if (c == '=') c = bam_is_rev(p->b) ? ',' : '.';
            else c = bam_is_rev(p->b) ? tolower(c) : toupper(c);
        }
//        putc(c, fp);
        tmpBase.push_back(c);
    }
//    else putc(p->is_refskip? (bam_is_rev(p->b)? '<' : '>') : '*', fp);
    if (p->indel > 0) {//insertion
//        putc('+', fp); printw(p->indel, fp);
        for (j = 1; j <= p->indel; ++j) {
            int c = seq_nt16_str[bam_seqi(bam_get_seq(p->b), p->qpos + j)];
//            putc(bam_is_rev(p->b)? tolower(c) : toupper(c), fp);
        }
    } else if (p->indel < 0) {//deletion
//        printw(p->indel, fp);
        for (j = 1; j <= -p->indel; ++j) {
            int c = (ref && (int) pos + j < ref_len) ? ref[pos + j] : 'N';
//            putc(bam_is_rev(p->b)? tolower(c) : toupper(c), fp);
        }
    }
//    if (p->is_tail) putc('$', fp);
}


#define MPLP_BCF        1
#define MPLP_VCF        (1<<1)
#define MPLP_NO_COMP    (1<<2)
#define MPLP_NO_ORPHAN  (1<<3)
#define MPLP_REALN      (1<<4)
#define MPLP_NO_INDEL   (1<<5)
#define MPLP_REDO_BAQ   (1<<6)
#define MPLP_ILLUMINA13 (1<<7)
#define MPLP_IGNORE_RG  (1<<8)
#define MPLP_PRINT_POS  (1<<9)
#define MPLP_PRINT_MAPQ (1<<10)
#define MPLP_PER_SAMPLE (1<<11)
#define MPLP_SMART_OVERLAPS (1<<12)

#ifdef __cplusplus
extern "C" {
#endif
void *bed_read(const char *fn);
void bed_destroy(void *_h);
int bed_overlap(const void *_h, const char *chr, int beg, int end);

#ifdef __cplusplus
}
#endif

typedef struct {
    char *ref[2];
    int ref_id[2];
    int ref_len[2];
} mplp_ref_t;

#define MPLP_REF_INIT {{NULL,NULL},{-1,-1},{0,0}}

typedef struct {
    samFile *fp;
    hts_itr_t *iter;
    bam_hdr_t *h;
    mplp_ref_t *ref;
    const mplp_conf_t *conf;
} mplp_aux_t;

typedef struct {
    int n;
    int *n_plp, *m_plp;
    bam_pileup1_t **plp;
} mplp_pileup_t;

static int mplp_get_ref(mplp_aux_t *ma, int tid, char **ref, int *ref_len) {
    mplp_ref_t *r = ma->ref;

    //printf("get ref %d {%d/%p, %d/%p}\n", tid, r->ref_id[0], r->ref[0], r->ref_id[1], r->ref[1]);

    if (!r || !ma->conf->fai) {
        *ref = NULL;
        return 0;
    }

    // Do we need to reference count this so multiple mplp_aux_t can
    // track which references are in use?
    // For now we just cache the last two. Sufficient?
    if (tid == r->ref_id[0]) {
        *ref = r->ref[0];
        *ref_len = r->ref_len[0];
        return 1;
    }
    if (tid == r->ref_id[1]) {
        // Last, swap over
        int tmp;
        tmp = r->ref_id[0];
        r->ref_id[0] = r->ref_id[1];
        r->ref_id[1] = tmp;
        tmp = r->ref_len[0];
        r->ref_len[0] = r->ref_len[1];
        r->ref_len[1] = tmp;

        char *tc;
        tc = r->ref[0];
        r->ref[0] = r->ref[1];
        r->ref[1] = tc;
        *ref = r->ref[0];
        *ref_len = r->ref_len[0];
        return 1;
    }

    // New, so migrate to old and load new
    free(r->ref[1]);
    r->ref[1] = r->ref[0];
    r->ref_id[1] = r->ref_id[0];
    r->ref_len[1] = r->ref_len[0];

    r->ref_id[0] = tid;
    r->ref[0] = faidx_fetch_seq(ma->conf->fai,
                                ma->h->target_name[r->ref_id[0]],
                                0,
                                INT_MAX,
                                &r->ref_len[0]);

    if (!r->ref[0]) {
        r->ref[0] = NULL;
        r->ref_id[0] = -1;
        r->ref_len[0] = 0;
        *ref = NULL;
        return 0;
    }

    *ref = r->ref[0];
    *ref_len = r->ref_len[0];
    return 1;
}

#ifdef __cplusplus
extern "C" {
#endif
extern int bam_realn(bam1_t *b, const char *ref);
extern int bam_prob_realn_core(bam1_t *b, const char *ref, int ref_len, int flag);
extern int bam_cap_mapQ(bam1_t *b, char *ref, int ref_len, int thres);

#ifdef __cplusplus
}
#endif

static int mplp_func(void *data, bam1_t *b) {

    char *ref;
    mplp_aux_t *ma = (mplp_aux_t *) data;
    int ret, skip = 0, ref_len;
    do {
        int has_ref;
        ret = ma->iter ? sam_itr_next(ma->fp, ma->iter, b) : sam_read1(ma->fp, ma->h, b);
        if (ret < 0) break;
        // The 'B' cigar operation is not part of the specification, considering as obsolete.
        //  bam_remove_B(b);
        if (b->core.tid < 0 || (b->core.flag & BAM_FUNMAP)) { // exclude unmapped reads
            skip = 1;
            continue;
        }
        if (ma->conf->rflag_require && !(ma->conf->rflag_require & b->core.flag)) {
            skip = 1;
            continue;
        }
        if (ma->conf->rflag_filter && ma->conf->rflag_filter & b->core.flag) {
            skip = 1;
            continue;
        }
        if (ma->conf->bed) { // test overlap
            skip = !bed_overlap(ma->conf->bed, ma->h->target_name[b->core.tid], b->core.pos, bam_endpos(b));
            if (skip) continue;
        }
        if (ma->conf->rghash) { // exclude read groups
            uint8_t *rg = bam_aux_get(b, "RG");
            skip = (rg && khash_str2int_get(ma->conf->rghash, (const char *) (rg + 1), NULL) == 0);
            if (skip) continue;
        }
        if (ma->conf->flag & MPLP_ILLUMINA13) {
            int i;
            uint8_t *qual = bam_get_qual(b);
            for (i = 0; i < b->core.l_qseq; ++i)
                qual[i] = qual[i] > 31 ? qual[i] - 31 : 0;
        }

        if (ma->conf->fai && b->core.tid >= 0) {
            has_ref = mplp_get_ref(ma, b->core.tid, &ref, &ref_len);
            if (has_ref && ref_len <= b->core.pos) { // exclude reads outside of the reference sequence
                fprintf(stderr, "[%s] Skipping because %d is outside of %d [ref:%d]\n",
                        __func__, b->core.pos, ref_len, b->core.tid);
                skip = 1;
                continue;
            }
        } else {
            has_ref = 0;
        }

        skip = 0;
        if (has_ref && (ma->conf->flag & MPLP_REALN))
            bam_prob_realn_core(b, ref, ref_len, (ma->conf->flag & MPLP_REDO_BAQ) ? 7 : 3);
        if (has_ref && ma->conf->capQ_thres > 10) {
            int q = bam_cap_mapQ(b, ref, ref_len, ma->conf->capQ_thres);
            if (q < 0) skip = 1;
            else if (b->core.qual > q) b->core.qual = q;
        }
        if (b->core.qual < ma->conf->min_mq) skip = 1;
        else if ((ma->conf->flag & MPLP_NO_ORPHAN) && (b->core.flag & BAM_FPAIRED) &&
                 !(b->core.flag & BAM_FPROPER_PAIR))
            skip = 1;
    } while (skip);
    return ret;
}


/*
 * Performs pileup
 * @param conf configuration for this pileup
 * @param n number of files specified in fn
 * @param fn filenames
 */
int SimplePileupViewer::SIMPLEmpileup(mplp_conf_t *conf, int n, char **fn) {
    bool BedEOF;
    hts_idx_t *idx = NULL;
    numBases = 0;
    avgDepth = 0;
    sdDepth = 0;

    mplp_aux_t **data;
    int i(0), tid(0), pos(0), *n_plp(0), beg0 = 0, end0 = INT_MAX, ref_len(0), max_depth(0), max_indel_depth(0);
    const bam_pileup1_t **plp;
    mplp_ref_t mp_ref = MPLP_REF_INIT;
    bam_mplp_t iter;
    bam_hdr_t *h = NULL; /* header of first file in input list */
    char *ref;
    void *rghash = NULL;
    FILE *pileup_fp = NULL;

    htsFile *bcf_fp = NULL;
    bcf_hdr_t *bcf_hdr = NULL;

    bam_sample_t *sm = NULL;
    kstring_t buf;
    mplp_pileup_t gplp;

    memset(&gplp, 0, sizeof(mplp_pileup_t));
    memset(&buf, 0, sizeof(kstring_t));
//    memset(&bc, 0, sizeof(bcf_call_t));
    data = (mplp_aux_t **) calloc(n, sizeof(mplp_aux_t *));
    plp = (const bam_pileup1_t **) calloc(n, sizeof(bam_pileup1_t *));
    n_plp = (int *) calloc(n, sizeof(int));
    sm = bam_smpl_init();

    if (n == 0) {
        fprintf(stderr, "[%s] no input file/data given\n", __func__);
        exit(EXIT_FAILURE);
    }

    // read the header of each file in the list and Initialize data
    for (i = 0; i < n; ++i) {
        bam_hdr_t *h_tmp;
        data[i] = (mplp_aux_t *) calloc(1, sizeof(mplp_aux_t));
        data[i]->fp = sam_open_format(fn[i], "rb", &conf->ga.in);
        if (!data[i]->fp) {
            fprintf(stderr, "[%s] failed to open %s: %s\n", __func__, fn[i], strerror(errno));
            exit(EXIT_FAILURE);
        }
        if (hts_set_opt(data[i]->fp, CRAM_OPT_DECODE_MD, 0)) {
            fprintf(stderr, "Failed to set CRAM_OPT_DECODE_MD value\n");
            exit(EXIT_FAILURE);
        }
        if (conf->fai_fname && hts_set_fai_filename(data[i]->fp, conf->fai_fname) != 0) {
            fprintf(stderr, "[%s] failed to process %s: %s\n",
                    __func__, conf->fai_fname, strerror(errno));
            exit(EXIT_FAILURE);
        }
        data[i]->conf = conf;
        data[i]->ref = &mp_ref;
        h_tmp = sam_hdr_read(data[i]->fp);
        if (!h_tmp) {
            fprintf(stderr, "[%s] fail to read the header of %s\n", __func__, fn[i]);
            exit(EXIT_FAILURE);
        }
        bam_smpl_add(sm, fn[i], (conf->flag & MPLP_IGNORE_RG) ? 0 : h_tmp->text);

        if (0) {
            SEQ_SM = std::string(sm->smpl[0]);//sm->n:# of sample; sm->m:max # of sample
            warning("The BAM or CRAM file contains more than 1 sample, please be aware the risk! Ignored...");
        } else {
            if (sm->n == 1) SEQ_SM = std::string(sm->smpl[0]);//sm->n:# of sample; sm->m:max # of sample
            else
                error("This BAM or CRAM file contains more than 1 sample, please demultiplex or separate first!");
        }
        idx = sam_index_load(data[i]->fp, fn[i]);
        if (idx == NULL) {
            fprintf(stderr, "[%s] fail to load index for %s\n", __func__, fn[i]);
            exit(EXIT_FAILURE);
        }
        region_t firstReg;
        conf->reg = new char[1024];
        REGET:
        firstReg = GetNextRegion(BedEOF);
        sprintf(conf->reg, "%s:%d-%d", firstReg.chr.c_str(), firstReg.beg + 1, firstReg.end);
        if (conf->reg) {
//            hts_idx_t *idx = sam_index_load(data[i]->fp, fn[i]);

            if ((data[i]->iter = sam_itr_querys(idx, h_tmp, conf->reg)) == 0) {
                fprintf(stderr, "[Warning::%s] initialization fail to parse region '%s' with %s, skip...\n", __func__,
                        conf->reg, fn[i]);
                if (BedEOF) {
                    delete[] conf->reg;
                    conf->reg = nullptr;
                    fprintf(stderr, "No reads found in any of the regions, exit!");
                    exit(EXIT_FAILURE);
                }
                goto REGET;
            }
            if (i == 0) beg0 = data[i]->iter->beg, end0 = data[i]->iter->end;
//            hts_idx_destroy(idx);
        } else
            data[i]->iter = NULL;

        if (i == 0) h = data[i]->h = h_tmp; // save the header of the first file
        else {
            // FIXME: check consistency between h and h_tmp
            bam_hdr_destroy(h_tmp);

            // we store only the first file's header; it's (alleged to be)
            // compatible with the i-th file's target_name lookup needs
            data[i]->h = h;
        }
    }
    // allocate data storage proportionate to number of samples being studied sm->n
    gplp.n = sm->n;
    gplp.n_plp = (int *) calloc(sm->n, sizeof(int));
    gplp.m_plp = (int *) calloc(sm->n, sizeof(int));
    gplp.plp = (bam_pileup1_t **) calloc(sm->n, sizeof(bam_pileup1_t *));

    fprintf(stderr, "[%s] %d sample(s) in %d input file(s)\n", __func__, sm->n, n);

    // init pileup
    iter = bam_mplp_init(n, mplp_func, (void **) data);
    if (conf->flag & MPLP_SMART_OVERLAPS) bam_mplp_init_overlaps(iter);
    max_depth = conf->max_depth;
    if (max_depth * sm->n > 1 << 20)
        fprintf(stderr, "[%s] Max depth is above 1M. Potential memory hog!\n", __func__);
    if (max_depth * sm->n < 8000) {
        max_depth = 8000 / sm->n;
        fprintf(stderr, "[%s] Set max per-file depth to %d\n", __func__, max_depth);
    }
    max_indel_depth = conf->max_indel_depth * sm->n;
    bam_mplp_set_maxcnt(iter, max_depth);
    bcf1_t *bcf_rec = bcf_init1();
    int ret;
    // begin pileup

    //adaptor for viewer
    int32_t globalIndex = 0;


    while (1) {
        ret = bam_mplp_auto(iter, &tid, &pos, n_plp, plp);
        //Fan changed:
        //Now we assume the bam files are sorted, then sites within each region should be consecutive,
        // then each time pos > end0, we can update beg0 and end0, as well as data[0]->iter

        //old implementation:
        //if (conf->reg && (pos < beg0 || pos >= end0)) continue; // out of the region requested
        //new implementation:
        if (ret <= 0 || pos >= end0) {
            if (BedEOF) break;
            char reg[1024];
            region_t tmp;
            REGET2:
            tmp = GetNextRegion(BedEOF);
            sprintf(reg, "%s:%d-%d", tmp.chr.c_str(), tmp.beg + 1, tmp.end);
            notice("Process %s...", reg);
            if ((data[0]->iter = sam_itr_querys(idx, data[0]->h, reg)) == 0) {
                fprintf(stderr, "[Warning::%s] fail to parse region '%s' with %s\n", __func__, reg, fn[0]);
                if (BedEOF) {
                    break;
                }
                goto REGET2;
            }
            beg0 = data[0]->iter->beg, end0 = data[0]->iter->end;
            bam_mplp_destroy(iter);
            iter = bam_mplp_init(n, mplp_func, (void **) data);
            if (conf->flag & MPLP_SMART_OVERLAPS) bam_mplp_init_overlaps(iter);
            continue;
        } else if (pos < beg0) continue;

/*if original -l option*/    if (conf->bed && tid >= 0 &&
                                 !bed_overlap(conf->bed, h->target_name[tid], pos, pos + 1))
            continue;

        mplp_get_ref(data[0], tid, &ref, &ref_len);
        //printf("tid=%d len=%d ref=%p/%s\n", tid, ref_len, ref, ref);
        {
            //fprintf(pileup_fp, "%s\t%d\t%c", h->target_name[tid], pos + 1, (ref && pos < ref_len)? ref[pos] : 'N');

            std::string chr = h->target_name[tid];
            bool existed(false);

            if (posIndex.find(chr) != posIndex.end())//chr existed
            {
                if (posIndex[chr].find(pos + 1) != posIndex[chr].end())//pos existed
                {
                    existed = true;
                } else {
                    posIndex[chr][pos + 1] = globalIndex;
                    globalIndex++;
                }
            } else {
                posIndex[chr][pos + 1] = globalIndex;
                globalIndex++;
            }

            if (existed) {
                warning("Duplicated marker %s:%d present, skip ...", chr.c_str(), pos + 1);
                continue;
            }

            std::vector<char> tmpBase, tmpQual;

            for (i = 0; i < n; ++i) {//for each bam file
                int j, cnt;
//                for (j = cnt = 0; j < n_plp[i]; ++j) {//each covered read
//                    const bam_pileup1_t *p = plp[i] + j;
//                    int c = p->qpos < p->b->core.l_qseq
//                            ? bam_get_qual(p->b)[p->qpos]
//                            : 0;
//                    if (c >= conf->min_baseQ) ++cnt;
//                }
                //fprintf(pileup_fp, "\t%d\t", cnt);
                if (n_plp[i] == 0) {// if no reads covered
                    //fputs("*\t*", pileup_fp);
                    //if (conf->flag & MPLP_PRINT_MAPQ) fputs("\t*", pileup_fp);
                    //if (conf->flag & MPLP_PRINT_POS) fputs("\t*", pileup_fp);
                    continue;
                } else {
                    /*calculate number of reads covering snps*/
//                    numBases += n_plp[i];
                    for (j = 0; j < n_plp[i]; ++j) {//each covered read in ith bam file
                        const bam_pileup1_t *p = plp[i] + j;
                        int c = p->qpos < p->b->core.l_qseq
                                ? bam_get_qual(p->b)[p->qpos]
                                : 0;
                        if (c >= conf->min_baseQ)//SimplePileupViewer Change
                        {
                            pileup_seq(pileup_fp, plp[i] + j, pos, ref_len, ref, tmpBase);
                        }
                    }
                    //putc('\t', pileup_fp);
                    for (j = 0; j < n_plp[i]; ++j) {
                        const bam_pileup1_t *p = plp[i] + j;
                        int c = p->qpos < p->b->core.l_qseq
                                ? bam_get_qual(p->b)[p->qpos]
                                : 0;
                        if (c >= conf->min_baseQ) {
                            c = c + 33 < 126 ? c + 33 : 126;
                            //putc(c, pileup_fp);
                            tmpQual.push_back(c);
                        }
                    }
//                    if (conf->flag & MPLP_PRINT_MAPQ) {//multiple pileups
//                        putc('\t', pileup_fp);
//                        for (j = 0; j < n_plp[i]; ++j) {
//                            const bam_pileup1_t *p = plp[i] + j;
//                            int c = bam_get_qual(p->b)[p->qpos];
//                            if ( c < conf->min_baseQ ) continue;
//                            c = plp[i][j].b->core.qual + 33;
//                            if (c > 126) c = 126;
//                            putc(c, pileup_fp);
//                        }
//                    }
//                    if (conf->flag & MPLP_PRINT_POS) {
//                        putc('\t', pileup_fp);
//                        int last = 0;
//                        for (j = 0; j < n_plp[i]; ++j) {
//                            const bam_pileup1_t *p = plp[i] + j;
//                            int c = bam_get_qual(p->b)[p->qpos];
//                            if ( c < conf->min_baseQ ) continue;
//
//                            if (last++) putc(',', pileup_fp);
//                            fprintf(pileup_fp, "%d", plp[i][j].qpos + 1); // FIXME: printf() is very slow...
//                        }
//                    }
                }
            }
            if(tmpBase.size()>0)
            {
                effectiveNumSite++;
                numBases +=tmpBase.size();
            }
            // putc('\n', pileup_fp);
            if (not existed) {
                baseInfo.push_back(tmpBase);
                qualInfo.push_back(tmpQual);
            }

        }

    }
    avgDepth = double(numBases) / GetNumMarker();
    hts_idx_destroy(idx);

    // clean up
//    free(bc.tmp.s);
    bcf_destroy1(bcf_rec);
    if (bcf_fp) {
        hts_close(bcf_fp);
        bcf_hdr_destroy(bcf_hdr);
//        bcf_call_destroy(bca);
//        free(bc.PL);
//        free(bc.DP4);
//        free(bc.ADR);
//        free(bc.ADF);
//        free(bc.fmt_arr);
//        free(bcr);
    }
    //if (pileup_fp && conf->output_fname) fclose(pileup_fp);
    bam_smpl_destroy(sm);
    free(buf.s);
    for (i = 0; i < gplp.n; ++i) free(gplp.plp[i]);
    free(gplp.plp);
    free(gplp.n_plp);
    free(gplp.m_plp);
//    bcf_call_del_rghash(rghash);
    bam_mplp_destroy(iter);
    bam_hdr_destroy(h);
    for (i = 0; i < n; ++i) {
        sam_close(data[i]->fp);
        if (data[i]->iter) hts_itr_destroy(data[i]->iter);
        free(data[i]);
    }
    free(data);
    free(plp);
    free(n_plp);
    free(mp_ref.ref[0]);
    free(mp_ref.ref[1]);
    return ret;
}

#define MAX_PATH_LEN 1024

int read_file_list(const char *file_list, int *n, char **argv[]) {
    char buf[MAX_PATH_LEN];
    int len, nfiles = 0;
    char **files = NULL;
    struct stat sb;

    *n = 0;
    *argv = NULL;

    FILE *fh = fopen(file_list, "r");
    if (!fh) {
        fprintf(stderr, "%s: %s\n", file_list, strerror(errno));
        return 1;
    }

    files = (char **) calloc(nfiles, sizeof(char *));
    nfiles = 0;
    while (fgets(buf, MAX_PATH_LEN, fh)) {
        // allow empty lines and trailing spaces
        len = strlen(buf);
        while (len > 0 && isspace(buf[len - 1])) len--;
        if (!len) continue;

        // check sanity of the file list
        buf[len] = 0;
        if (stat(buf, &sb) != 0) {
            // no such file, check if it is safe to print its name
            int i, safe_to_print = 1;
            for (i = 0; i < len; i++)
                if (!isprint(buf[i])) {
                    safe_to_print = 0;
                    break;
                }
            if (safe_to_print)
                fprintf(stderr, "The file list \"%s\" appears broken, could not locate: %s\n", file_list, buf);
            else
                fprintf(stderr, "Does the file \"%s\" really contain a list of files and do all exist?\n", file_list);
            return 1;
        }

        nfiles++;
        files = (char **) realloc(files, nfiles * sizeof(char *));
        files[nfiles - 1] = strdup(buf);
    }
    fclose(fh);
    if (!nfiles) {
        fprintf(stderr, "No files read from %s\n", file_list);
        return 1;
    }
    *argv = files;
    *n = nfiles;
    return 0;
}

#undef MAX_PATH_LEN

#include <iostream>
#include <fstream>
#include <sstream>

SimplePileupViewer::SimplePileupViewer() {

}

SimplePileupViewer::SimplePileupViewer(std::vector<region_t> *BedPtr, const char *bamFile, const char *faiFile,
                                       const char *bedFile, int nfiles) {
    bedVec = BedPtr;
    regIndex = 0;
    int c;
    const char *file_list = bamFile;
    char **fn = NULL;
//    int nfiles = 0, use_orphan = 0;
//    mplp_conf_t mplp;
    memset(&mplp, 0, sizeof(mplp_conf_t));
    mplp.min_mq = 13;
    mplp.min_baseQ = 2;
    mplp.capQ_thres = 40;
    mplp.max_depth = 250;
    mplp.max_indel_depth = 250;
    mplp.openQ = 40;
    mplp.extQ = 20;
    mplp.tandemQ = 100;
    mplp.min_frac = 0.002;
    mplp.min_support = 1;
    mplp.flag = /*MPLP_NO_ORPHAN |*/ MPLP_REALN | MPLP_SMART_OVERLAPS;
    mplp.argc = 0;
    mplp.argv = 0;
    mplp.rflag_filter = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP;
    mplp.output_fname = NULL;
    sam_global_args_init(&mplp.ga);
    if (faiFile == NULL) {
        std::cerr << "fai file open failed!" << std::endl;
        exit(EXIT_FAILURE);
    } else {
        mplp.fai = fai_load(faiFile);
        mplp.fai_fname = const_cast<char *>(faiFile);
        if (mplp.fai == NULL) {
            std::cerr << "open fai file failed!";
            exit(EXIT_FAILURE);
        }
        mplp.ga.reference = mplp.fai_fname;
    }

    // In the original version the whole BAM was streamed which is inefficient
    //  with few BED intervals and big BAMs. Todo: devise a heuristic to determine
    //  best strategy, that is streaming or jumping.
    mplp.bed = bed_read(bedFile);
    if (!mplp.bed) {
        print_error_errno("mpileup", "Could not read file \"%s\"", optarg);
        exit(EXIT_FAILURE);
    }
    int ret(0);
    if (nfiles > 1) {
        if (read_file_list(file_list, &nfiles, &fn)) {
            std::cerr << "open bam file list failed!" << std::endl;
            exit(EXIT_FAILURE);
        }
        ret = SIMPLEmpileup(&mplp, nfiles, fn);
        for (c = 0; c < nfiles; c++) free(fn[c]);
        free(fn);
    } else {
        fn = new char *[1];
        fn[0] = strdup(bamFile);
        ret = SIMPLEmpileup(&mplp, nfiles, fn);
        free(fn[0]);
        delete[] fn;
    }
    if (mplp.fai) fai_destroy(mplp.fai);
    if (mplp.bed) bed_destroy(mplp.bed);
}

SimplePileupViewer::SimplePileupViewer(const BED& BedFromEstimator, const std::string &pileupFile):bedTable(BedFromEstimator){
    ReadPileup(pileupFile);
}

static std::string ParsePileupSeq(std::string seq, std::string refAllele)
{
    std::string newSeq;
    std::transform(seq.begin(), seq.end(), seq.begin(),
                   [](unsigned char c){ return std::toupper(c); }
    );
    for(int i=0; i!=seq.size(); i++)
    {
        if(seq[i]=='+' or seq[i]=='-')
        {
            int tmpIndex=i+1;
            while(tmpIndex != seq.size() and std::isdigit(seq[tmpIndex]))
            {
                tmpIndex++;
            }
            int digitLen=tmpIndex-(i+1);
            int clipLen=std::stoi(seq.substr(i+1,digitLen));
            i+= digitLen + clipLen;
        }
        else if(seq[i] == '^')
        {
            i+=1;
        }
        else if(seq[i]=='.' or seq[i]==',') newSeq += refAllele[0];
        else if(seq[i] == 'A' or seq[i]== 'G' or seq[i] =='C' or seq[i] == 'T' or seq[i]=='N') newSeq+=seq[i];
    }
    return newSeq;
}
int SimplePileupViewer::ReadPileup(const std::string &filePath) {

    int globalIndex=0;

    //pileup variables
    std::string pChr;
    int pPos;
    std::string refAllele;
    int depth;
    std::string seq;
    std::string qual;
    std::ifstream fin(filePath);
    numBases = 0;
    if(not fin.is_open())
    {
        std::cerr<<"open file "<<filePath<<" failed!"<<std::endl;
        exit(EXIT_FAILURE);
    }
    std::string pileupLine;
    while(getline(fin,pileupLine)) {

        std::stringstream ss(pileupLine);
        ss>>pChr>>pPos>>refAllele>>depth>>seq>>qual;
        if(seq.find_first_of(".,")!=std::string::npos)//using ., format
        {
            if(refAllele==".")
            {
                std::cerr<<"Pileup format error: cannot find ref allele, exit!\n";
                exit(EXIT_FAILURE);
            }
        }
//        seq=ParsePileupSeq(seq,refAllele);

        if(bedTable.find(pChr)==bedTable.end())
            continue;
        else if(bedTable[pChr].find(pPos)==bedTable[pChr].end())
            continue;


        int tmpIndex(0);
        bool existed(false);

        if (posIndex.find(pChr) != posIndex.end())//chr existed
        {
            if (posIndex[pChr].find(pPos) != posIndex[pChr].end())//pos existed
            {
                tmpIndex = posIndex[pChr][pPos];
                existed = true;
            } else {
                posIndex[pChr][pPos] = globalIndex;
                globalIndex++;
            }
        } else {
            posIndex[pChr][pPos] = globalIndex;
            globalIndex++;
        }

        std::vector<char> tmpBase, tmpQual;
        if (existed) {
            std::cerr<<"[WARNING] The pileup file has duplicated lines! Merged here"<<std::endl;
            tmpBase = GetBaseInfoAt(pChr, pPos);
            tmpQual = GetQualInfoAt(pChr, pPos);
        }

        std::copy( seq.begin(), seq.end(), std::back_inserter(tmpBase));
        std::copy( qual.begin(), qual.end(), std::back_inserter(tmpQual));

        if (not existed) {
            baseInfo.push_back(tmpBase);
            qualInfo.push_back(tmpQual);
        }
        numBases += depth;
        depth = 0;
        seq = "";
        qual = "";
        effectiveNumSite++;
    }
    avgDepth = (double)numBases/GetNumMarker();
    return 0;
}

SimplePileupViewer::~SimplePileupViewer() {

}
