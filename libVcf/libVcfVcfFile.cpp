#include <math.h>
#include <limits.h>

#include "libVcfVcfFile.h"
using namespace libVcf;

std::vector<double> VcfHelper::vPhred2Err;
StringArray VcfHelper::asChromNames;
std::vector<int> VcfHelper::vnChromNums;
int VcfHelper::char2twobit[255];

bool dummy1 = VcfHelper::initPhred2Error();
bool dummy2 = VcfHelper::initChromNamesNums();
bool dummy3 = VcfHelper::initChar2TwoBit();

bool VcfHelper::initPhred2Error(int maxPhred) {
  vPhred2Err.clear();
  for(int i=0; i <= maxPhred; ++i) {
    vPhred2Err.push_back(pow(0.1, i/10.));
  }
  return true;
}

bool VcfHelper::initChromNamesNums() {
  asChromNames.Add("X");
  vnChromNums.push_back(23);
  asChromNames.Add("Y");
  vnChromNums.push_back(24);
  asChromNames.Add("XY");
  vnChromNums.push_back(25);
  asChromNames.Add("MT");
  vnChromNums.push_back(26);
  return true;
}

bool VcfHelper::initChar2TwoBit() {
  memset(char2twobit, sizeof(char), 255);
  char2twobit['A'] = 0;
  char2twobit['C'] = 1;
  char2twobit['G'] = 2;
  char2twobit['T'] = 3;
  return true;
}

uint64_t VcfHelper::str2TwoBits(const char* s, int len) {
  if ( len > 32 ) {
    throw VcfFileException("Cannot encode sequences of length > %d",len);
  }
  uint64_t ret = 0;
  for(int i=0; i < len; ++i) {
    ret = (ret << 2); // shift two bits
    ret |= ( char2twobit[(int)s[i]] & 0x03 ); // and take ORs
  }
  return ret;
}

int VcfHelper::chromName2Num(const String & chr) {
  int n = atoi(chr.c_str());
  if ( n > 0 ) { return n; }
  else {
    String s = chr;
    if ( s.Left(3).Compare("chr") == 0 ) {
      n = atoi(s.SubStr(3).c_str());
      if ( n > 0 ) { return n; }
      s = s.SubStr(3);
    }
    for(int i=0; i < asChromNames.Length(); ++i) {
      if ( s.Compare(asChromNames[i]) == 0 ) {
	return vnChromNums[i];
      }
    }
  }
  throw VcfFileException("Cannot recognize chromosome %s",chr.c_str());
}

int VcfHelper::compareGenomicPos(const String & chr1, int pos1, const String& chr2, int pos2) {
  //Logger::gLogger->writeLog("VcfHelper::compareGenomicPos %s:%d - %s:%d",chr1.c_str(), pos1, chr2.c_str(), pos2);
  int nchr1 = chromName2Num(chr1);
  int nchr2 = chromName2Num(chr2);

  if ( nchr1 == nchr2 ) {
    return pos1 - pos2;
  }
  else {
    return (nchr1 - nchr2)*1000000;
  }
}

VcfInd::VcfInd(const String& indID, const String& famID, const String& fatID, const String& motID, const String& gender) {
  sIndID = indID;
  sFamID = famID;
  sFatID = (fatID.Length() == 0) ? "0" : fatID;
  sMotID = (motID.Length() == 0) ? "0" : motID;
  if ( gender.Compare("1") == 0 ) {
    this->gender = MALE;
  }
  else if ( gender.Compare("2") == 0 ) {
    this->gender = FEMALE;
  }
  else {
    this->gender = UNKNOWN;
  }
}

VcfFile::VcfFile() {
  iFile = NULL;
  nBuffers = 0;
  nNumMarkers = 0;
  nHead = 0;
  bSiteOnly = false;
  bParseGenotypes = true;
  bParseDosages = true;
  bParseValues = true;
  bUpgrade = false;
  bEOF = false;
  nMinGD = 0;
  nMinGQ = 0;
}

VcfFile::~VcfFile() {
  reset();
}

void VcfFile::reset() {
  if ( iFile != NULL ) 
    ifclose(iFile);
  iFile = NULL;
  for(int i=0; i < (int) vpVcfInds.size(); ++i) {
    delete vpVcfInds[i];
  }
  vpVcfInds.clear();
  for(int i=0; i < (int) vpVcfMarkers.size(); ++i) {
    delete vpVcfMarkers[i];
  }
  vpVcfMarkers.clear();
  asMetaKeys.Clear();
  asMetaValues.Clear();

  nNumLines = 0;
}

void VcfFile::openForRead(const char* filename, int nbuf) {
  reset();
  
  iFile = ifopen(filename,"rb");
  if ( iFile == NULL ) {
    throw VcfFileException("Failed opening file %s - %s",filename, strerror(errno));
  }
  nBuffers = nbuf;
  nNumMarkers = 0;
  nHead = 0;
  if ( nBuffers == 0 ) { // infinite buffer size
    // do not set size of markers
  }
  else {
    vpVcfMarkers.resize( nBuffers );
    for(int i=0; i < nBuffers; ++i) {
      VcfMarker* p = new VcfMarker;
      vpVcfMarkers[i] = p;
    }
  }
  parseMeta();
  parseHeader();

  if ( bUpgrade ) {
    upgradeMetaLines();
  }
}

int VcfFile::readLine() {
  int retval;
  retval = line.ReadLine(iFile);
  if ( retval > 0 ) ++nNumLines;
  return retval;
}

// copied from Matthew Flickenger/Snyder
void VcfFile::parseMeta() {
  do {
    if ( readLine() <= 0 ) break;
    if ( ifeof(iFile) || !isMetaLine() ) break;
    
    parseMetaLine();
  } while (1);
  verifyMetaLines();
}

void VcfFile::parseMetaLine() {
  int equalPos = line.FindChar('=');
  if ( equalPos > 0 ) {
    asMetaKeys.Add( line.Mid(2,equalPos-1) );
    asMetaValues.Add( line.SubStr(equalPos+1) );
  }
  else {
    asMetaKeys.Add( line.SubStr(2) );
    asMetaValues.Add ("");
  }
}

// copied from Matthew Flickenger/Snyder
void VcfFile::verifyMetaLines() {
  if ( asMetaKeys.Length() < 1 ) {
    //throw VcfFileException("No meta ## lines found. ##format is required");
  }
  if ( asMetaKeys.Find("format") == -1 && asMetaKeys.Find("fileformat") == -1 ) {
    //throw VcfFileException("Required ##format line not found");
  }
}

// specific to upgrade from VCF 3.3 to VCF 4.0
void VcfFile::upgradeMetaLines() {
  std::vector<int> toDelete;

  for(int i=0; i < asMetaKeys.Length(); ++i) {
    if ( asMetaKeys[i].Compare("fileformat") == 0 ) {
      if ( asMetaValues[i].Compare("VCFv3.3") == 0 ) {
	asMetaValues[i] = "VCFv4.0";
      }
    }
    else if ( asMetaKeys[i].Compare("FILTER") == 0 ) {
      toDelete.push_back(i);
    }
    else if ( asMetaKeys[i].Compare("FORMAT") == 0 ) {
      StringArray tok;
      tok.ReplaceColumns(asMetaValues[i],',');
      asMetaValues[i].printf("<ID=%s,Number=%s,Type=%s,Description=%s>",tok[0].c_str(),tok[2].c_str(),tok[1].c_str(),tok[3].c_str());
    }
  }

  for(int i=0; i < (int)toDelete.size(); ++i) {
    asMetaValues.Delete(toDelete[toDelete.size()-i-1]);
    asMetaKeys.Delete(toDelete[toDelete.size()-i-1]);
  }

  // add INFO and FILTER?
}

void VcfFile::parseHeader() {
  if ( isHeaderLine() ) {
    lineTokens.ReplaceColumns(line, '\t');
    verifyHeaderLine();
    if ( lineTokens.Length() > 8 ) {
      int sampleCount = lineTokens.Length() - 9;
      vpVcfInds.resize(sampleCount);
      for(int i=0; i < sampleCount; ++i) {
	vpVcfInds[i] = new VcfInd(lineTokens[9+i]);
      }
    }
    else {
      for(int i=0; i < (int) vpVcfInds.size(); ++i) {
	delete vpVcfInds[i];
      }
      vpVcfInds.clear();
    }
  }
  else {
    throw VcfFileException("Header line is not found : #CHROM...");
  }
}

void VcfFile::verifyHeaderLine() {
   String expected[9] = {"#CHROM", "POS","ID","REF",
      "ALT","QUAL","FILTER","INFO", "FORMAT"};
   for(int i=0; i<8; i++) {
      if(lineTokens[i] != expected[i]) {
	throw VcfFileException("Bad header: expected %s, found %s",expected[i].c_str(),lineTokens[i].c_str());
      }
   }
   if (lineTokens.Length() > 9) {
      if(lineTokens[8] != expected[8]) {
	throw VcfFileException("Bad header: expected %s, found %s, token Length == %d",expected[8].c_str(),lineTokens[8].c_str(),lineTokens.Length());
      }
   }
}

bool VcfFile::isMetaLine() {
   return line[0] == '#' && line[1] == '#';
}

bool VcfFile::isHeaderLine() {
   return line[0] == '#' && line[1] == 'C';
}

int VcfFile::getMetaCount() {
   return asMetaKeys.Length();
}

int VcfFile::getMetaIndex(const String& meta) {
   return asMetaKeys.Find(meta);
}

bool VcfFile::hasMetaKey(const String& meta) {
   return asMetaKeys.Find(meta)!=-1;
}

String VcfFile::getMetaKey(int offset) {
   return asMetaKeys[offset];
}

String VcfFile::getMetaValue(int offset) {
   return asMetaValues[offset];
}
String VcfFile::getMetaValue(const String& meta) {
   return asMetaValues[getMetaIndex(meta)];
}

String VcfFile::getMetaValue(int offset, const String& ifmissing) {
   if (offset > -1) {
      return getMetaValue(offset);
   } else {
      return ifmissing;
   }
}
String VcfFile::getMetaValue(const String& meta, const String& ifmissing) {
   return getMetaValue(getMetaIndex(meta), ifmissing);
}

int VcfFile::getSampleCount() {
  return static_cast<int>(vpVcfInds.size());
}

int VcfFile::getSampleIndex(const String& name) {
  for(int i=0; i < (int)vpVcfInds.size(); ++i) {
    if ( name.Compare(vpVcfInds[i]->sIndID) == 0 ) {
      return i;
    }
  }
  return -1;
}

bool VcfFile::hasSample(String name) {
   return getSampleIndex(name) != -1;
}

String VcfFile::getSampleID(int offset) {
  if ( offset < (int)vpVcfInds.size() ) {
    return vpVcfInds[offset]->sIndID;
  }
  else {
    throw VcfFileException("VcfFile::getSampleID(%d) - Array out of bound (>%d)",offset,(int)vpVcfInds.size());
  }
}

BedFile::BedFile() {
  iFile = NULL;
  iBimFile = NULL;
  iFamFile = NULL;
  bAllowFlip = true; // allow flip
  bRefIsAllele1 = true;
  pBedBuffer = NULL;
  nBytes = 0;
}

BedFile::~BedFile() {
  if ( iBimFile != NULL ) {
    ifclose(iBimFile);
  }
  iBimFile = NULL;

  if ( iFamFile != NULL ) {
    ifclose(iFamFile);
  }
  iFamFile = NULL;
  if ( iFile != NULL ) {
    ifclose(iFile);
  }
  if ( pBedBuffer != NULL ) {
    delete[] pBedBuffer;
  }
}

void BedFile::openForRead(const char* bfile, const char* reffile, int nbuf) {
  String s = String(bfile);
  String bedFile = s + ".bed";
  String bimFile = s + ".bim";
  String famFile = s + ".fam";
  openForRead(bedFile.c_str(), bimFile.c_str(), famFile.c_str(), reffile, nbuf);
}

void BedFile::openForRead(const char* bedFile, const char* bimFile, const char* famFile, const char* refFile, int nbuf) {
  StringArray tokens;

  reset();

  iFile = ifopen(bedFile,"rb");
  if ( iFile == NULL ) {
    throw VcfFileException("Failed opening file %s - %s",bedFile,strerror(errno));
  }
  
  // read magic numbers
  char magicNumbers[3] = {0x6c,0x1b,0x01};
  char firstThreeBytes[3];
  ifread( iFile, firstThreeBytes, 3 );
  for(int i=0; i < 3; ++i) {
    if ( firstThreeBytes[i] != magicNumbers[i] ) {
      throw VcfFileException("The magic numbers do not match in BED file %s",bedFile);
    }
  }

  iBimFile = ifopen(bimFile,"rb");
  iFamFile = ifopen(famFile,"rb");
  sRefFile = refFile;

  while( 1 ) {
    int ret = line.ReadLine(iFamFile);
    if ( ret <= 0 ) break;
    tokens.ReplaceTokens(line, " \t\r\n");
    if ( tokens.Length() < 5 ) {
      throw VcfFileException("Less then 5 columns are observed in FAM file");
    }
    VcfInd* p = new VcfInd(tokens[1],tokens[0],tokens[2],tokens[3],tokens[4]);
    vpVcfInds.push_back(p);
  }

  //Logger::gLogger->writeLog("Finished loading %d individuals from FAM file",(int)vpVcfInds.size());

  nBytes = (vpVcfInds.size()+3)/4;
  if ( pBedBuffer != NULL ) { delete[] pBedBuffer; }
  pBedBuffer = new char[nBytes];

  nBuffers = nbuf;
  nNumMarkers = 0;
  nHead = 0;

  bParseGenotypes = true;
  bParseDosages = false;
  bParseValues = false;

  if ( nBuffers == 0 ) { // infinite buffer size
    // do not set size of markers
  }
  else {
    vpVcfMarkers.resize( nBuffers );
    for(int i=0; i < nBuffers; ++i) {
      VcfMarker* p = new VcfMarker;
      vpVcfMarkers[i] = p;
    }
  }

  genomeSequence.setReferenceName(sRefFile.c_str());
  genomeSequence.useMemoryMap(true);

  //Logger::gLogger->writeLog("Loading reference file %s",sRefFile.c_str());

  if ( genomeSequence.open() ) {
    // write a message that new index file is being created
    if ( genomeSequence.create(false) ) {
      throw VcfFileException("Failed creating index file of the reference. Please check the file permission");
    }
    if ( genomeSequence.open() ) {
      throw VcfFileException("Failed opening index file of the reference.");
    }
  }

  // set a default MetaLines
  asMetaKeys.Add("fileformat");
  asMetaValues.Add("VCFv4.0");
  asMetaKeys.Add("INFO");
  asMetaValues.Add("<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">");
  asMetaKeys.Add("INFO");
  asMetaValues.Add("<ID=AC,Number=.,Type=Integer,Description=\"Allele Count\">");
  asMetaKeys.Add("INFO");
  asMetaValues.Add("<ID=AN,Number=1,Type=Integer,Description=\"Number of Alleles With Data\">");
  asMetaKeys.Add("INFO");
  asMetaValues.Add("<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency\">");
  asMetaKeys.Add("FORMAT");
  asMetaValues.Add("<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
  asMetaKeys.Add("FORMAT");
  asMetaValues.Add("<ID=PL,Number=.,Type=Integer,Description=\"Phred-scale Genotype Likelihood\">");
}


bool VcfFile::iterateMarker() {
  if ( readLine() <= 0 ) { 
    bEOF = true;
    return false; 
  }

  lineTokens.ReplaceColumns(line, '\t');

  VcfMarker* pMarker;

  if ( nBuffers == 0 ) {
    pMarker = new VcfMarker;
    vpVcfMarkers.push_back(pMarker);
    ++nNumMarkers;
    ++nHead;
  }
  else {
    // make a circular list with constant size nBuffer
    pMarker = vpVcfMarkers[nHead];
    if ( pMarker == NULL ) {
      throw VcfFileException("nHead = %d, nBuffers = %d, pMarker is NULL",nHead,nBuffers);
    }
    nHead = (nHead+1) % nBuffers;
    ++nNumMarkers;
  }

  //Logger::gLogger->writeLog("foo, %s:%s",lineTokens[0].c_str(),lineTokens[1].c_str());

  try {
    pMarker->setChrom(lineTokens[0]);
    pMarker->setPos(lineTokens[1]);
    pMarker->setID(lineTokens[2]);
    pMarker->setRef(lineTokens[3]);
    pMarker->setAlts(lineTokens[4]);
    pMarker->setQual(lineTokens[5]);
    pMarker->setFilters(lineTokens[6]);

    if ( ( lineTokens.Length() >= 9 )  && ( !bSiteOnly ) ) {
      pMarker->setFormat(lineTokens[8], bUpgrade);

      int offset = 9;
      if ( lineTokens[9].IsEmpty() ) { // For handling bug in glfMultiples
	++offset;
      }

      pMarker->setSampleSize(lineTokens.Length()-offset, bParseGenotypes, bParseDosages, bParseValues);
      for(int i=offset; i < lineTokens.Length(); ++i) {
	pMarker->setSample(i-offset, lineTokens[i], bParseGenotypes, bParseDosages, bParseValues, nMinGD, nMinGQ);
      }
    }
    pMarker->setInfo(lineTokens[7], bUpgrade, ( ( nMinGD > 1 ) || ( nMinGQ > 0 ) ) );
  }
  catch (VcfFileException exc) {
    // add the line number to the error message
    throw VcfFileException(exc.msg + " See line " + nNumLines + ".");
  }
  return true;
}

void VcfMarker::setChrom(const String& s) {
  sChrom = s;
}

void VcfMarker::setPos(const String& s) {
  nPos = atoi(s.c_str());
}

void VcfMarker::setID(const String& s) {
  sID = s;
}

void VcfMarker::setRef(const String& s) {
  sRef = s;
  sRef.ToUpper();
}

void VcfMarker::setAlts(const String& s) {
  asAlts.Clear();
  if ( s.FindChar(',') == -1 ) {
    asAlts.Add(s);
    asAlts[0].ToUpper();
  }
  else {
    asAlts.ReplaceColumns(s,',');
    for(int i=0; i < asAlts.Length(); ++i) {
      asAlts[i].ToUpper();
    }
  }
}

void VcfMarker::setQual(const String& s) {
  if ( s.Compare(".") == 0 ) {
    fQual = -1;
  }
  else {
    fQual = atof(s.c_str());
  }
}

void VcfMarker::setFilters(const String& s) {
  asFilters.ReplaceColumns(s,';');
}

void VcfMarker::setInfo(const String& s, bool upgrade, bool updateAC) {
  if ( s[0] == '.' ) {
    if ( asInfoKeys.Length() > 0 ) {
      asInfoKeys.Clear();
      asInfoValues.Clear();
      bPreserved = false;
    }
    return;
  }
  //fprintf(stderr,"looking new line:\n");
  tmpTokens.ReplaceColumns(s, ';');

  if ( tmpTokens.Length() != asInfoKeys.Length() ) {
    bPreserved = false;

    asInfoKeys.Clear();
    asInfoValues.Clear();
    asInfoKeys.Dimension(tmpTokens.Length());
    asInfoValues.Dimension(tmpTokens.Length());

    for(int i=0; i < tmpTokens.Length(); ++i) {
      int equalsPos = tmpTokens[i].FindChar('=');
      if ( equalsPos == -1 ) {
	asInfoKeys[i] = tmpTokens[i];
	asInfoValues[i] = "";
      }
      else {
	asInfoKeys[i] = tmpTokens[i].Mid(0,equalsPos-1);
	asInfoValues[i] = tmpTokens[i].Mid(equalsPos+1,tmpTokens[i].Length()-1);
      }
    }
  }
  else {
    for(int i=0; i < tmpTokens.Length(); ++i) {
      int equalsPos = tmpTokens[i].FindChar('=');
      if ( equalsPos == -1 ) {
	if ( tmpTokens[i].Compare(asInfoKeys[i]) != 0 ) {
	  bPreserved = false;
	  asInfoKeys[i] = tmpTokens[i];
	}
	asInfoValues[i] = "";
      }
      else {
		 // fprintf(stderr, "compare %s:%s\t", tmpTokens[i].c_str(), asInfoKeys[i].c_str());
	if ( strncmp(tmpTokens[i].c_str(),asInfoKeys[i].c_str(),equalsPos) != 0 )  {// BUG!!!!!!!equalsPos should not minus 1
	  bPreserved = false;
	  asInfoKeys[i] = tmpTokens[i].Mid(0,equalsPos-1);
	}
	asInfoValues[i] = tmpTokens[i].Mid(equalsPos+1,tmpTokens[i].Length()-1);      
      }
    }
  }
  //fprintf(stderr, "end looking new line:\n");
  // upgrade INFO field entries from glfMultiples 06/16/2010 (VCFv3.3) to VCFv4.0 format
  if ( upgrade ) {
    if ( asInfoKeys[0].Compare("depth") == 0 ) {
      asInfoKeys[0] = "DP";
    }
    if ( asInfoKeys[1].Compare("mapQ") == 0 ) {
      asInfoKeys[1] = "MQ";
    }

    bool bMAF = false;
    double AF = atof(asInfoValues[2].c_str());
    if ( asInfoKeys[2].Compare("MAF") == 0 ) {
      bMAF = true;
      asInfoKeys[2] = "AF";
    }
    else if ( asInfoKeys[2].Compare("AF") == 0 ) {
      bMAF = false;
    }
    else {
      throw VcfFileException("VcfMarker::setInfo() : Cannot upgrade info field entry %s",s.c_str());
    }

    // calculate NS, AC, AN, AB statistics, assuming that sampleValues are already set
    int n = getSampleSize();
    int m = asFormatKeys.Length();
    //StringArray asPLs;
    if ( ( n * m == asSampleValues.Length() ) && ( m == 4 ) && ( n > 0 ) ) {
      int NS = 0;
      int ACs[3] = {0,0,0};
      std::vector<int> vnPLs;
      std::vector<int> vnDPs;
      vnPLs.resize(n * 3);
      vnDPs.resize(n);

      bool multiAlt = (asAlts.Length() > 1) ? true : false;

      for(int i=0; i < n; ++i) {
	int g1 = asSampleValues[i*m][0] - '0';
	int g2 = asSampleValues[i*m][2] - '0';
	int dp = atoi(asSampleValues[i*m+1].c_str());
	vnDPs[i] = dp;

	const char* s = asSampleValues[i*m+3].c_str();

	// read the GL fields to assign vnPL values for calculating allele balanace
	int pl = 0;
	int k = 0;
	for(int j=0; s[j] != '\0'; ++j) {
	  if ( s[j] == ',' ) {
	    if ( multiAlt ) {
	      if ( k == 2 ) { vnPLs[3*i+0] = pl; }
	      else if ( k == 4 ) { vnPLs[3*i+1] = pl; }
	      //else { fprintf(stderr,"%s\n",s); abort(); }
	    }
	    else {
	      vnPLs[3*i+k] = pl;
	    }
	    ++k;
	    pl = 0;
	  }
	  else {
	    pl = pl * 10 + (s[j]-'0');
	  }
	}
	if ( multiAlt ) {
	  if ( k == 5 ) { vnPLs[3*i+2] = pl; }
	  else { fprintf(stderr,"%s, k=%d\n",s,k); abort(); }
	}
	else {
	  vnPLs[3*i+k] = pl;
	}	

	//for(int j=0, k=3*i; s[j] != '\0'; ++j) {
	//  if ( s[j] == ',' ) {
        //    ++k;
	//  }
	//  else {
	//    vnPLs[k] = vnPLs[k] * 10 + (s[j]-'0');
	//  }
	//}

	/*
	asPLs.ReplaceColumns(asSampleValues[i*m+3],',');
	vnPLs[3*i+0] = atoi(asPLs[0].c_str());
	vnPLs[3*i+1] = atoi(asPLs[1].c_str());
	vnPLs[3*i+2] = atoi(asPLs[2].c_str());
	*/

	if ( dp > 0 ) {
	  ++NS;
	  ++ACs[g1];
	  ++ACs[g2];
	}
      }

      if ( bMAF ) {  // if bMAF is true, flip AF base on AC/NS
	if ( ACs[2] > 0 ) { // triallelic 
	  if ( NS < ACs[2] ) { // AF is greater than 0.5
	    AF = 1. - AF;
	  }
	}
	else { // biallelic 
	  if ( NS < ACs[1] ) {
	    AF = 1. - AF;
	  }
	}
      }

      if ( bMAF ) {
	asInfoValues[2].printf("%.3lf",AF);
      }

      // bound allele frequency estimates
      if ( AF < 1e-6 ) {
	AF = 1e-6;
      }
      else if ( 1.-AF < 1e-6 ) {
	AF = 1.-1e-6;
      }

      // from here Tom Blackwell's routine, documented at
      // http://genome.sph.umich.edu/wiki/Genotype_Likelihood_Based_Allele_Balance
      double ABnum = 0.5;
      double ABden = 1.0;

      double GPs[3];
      GPs[0] = (1.-AF)*(1.-AF);
      GPs[1] = 2.*AF*(1.-AF);
      GPs[2] = AF*AF;

      for(int i=0; i < n; ++i) {
	int nrefDen = vnPLs[i*3+2]+vnPLs[i*3+0]-2*vnPLs[3*i+1]+6*vnDPs[i];
	if ( nrefDen < 4 ) {
	  nrefDen = 4;
	}
	if ( nrefDen < abs(vnPLs[3*i+0]-vnPLs[3*i+2]) ) {
	  nrefDen = abs(vnPLs[3*i]-vnPLs[3*i+2]);
	}

	double nref = 0.5 * vnDPs[i] * (1.0 + (double)(vnPLs[i*3+2]-vnPLs[i*3+0])/(double)nrefDen);
	double pHet = VcfHelper::vPhred2Err[vnPLs[3*i+1]]*GPs[1]/(GPs[0]*VcfHelper::vPhred2Err[vnPLs[3*i]] + GPs[1]*VcfHelper::vPhred2Err[vnPLs[3*i+1]] + GPs[2]*VcfHelper::vPhred2Err[vnPLs[3*i+2]]);
	ABnum += (pHet * nref);
	ABden += (pHet * vnDPs[i]);
      }
      double AB = ABnum/ABden;

      String sNS, sAC, sAN, sAB;
      sNS.printf("%d",NS);
      if ( ACs[2] > 0 ) {
	sAC.printf("%d,%d",ACs[1],ACs[2]);
	asInfoValues[2].printf("%.3lf,%.3lf",1.-AF,AF);
      }
      else {
	sAC.printf("%d",ACs[1]);
      }

      sAN.printf("%d",2*NS);
      sAB.printf("%.4lf",AB);

      asInfoKeys.Add("NS");
      asInfoKeys.Add("AC");
      asInfoKeys.Add("AN");
      asInfoKeys.Add("AB");
      asInfoValues.Add(sNS);
      asInfoValues.Add(sAC);
      asInfoValues.Add(sAN);
      asInfoValues.Add(sAB);
    }
  }
  else if ( updateAC ) {
    int ACindex = -1;
    int ANindex = -1;
    int NSindex = -1;

    ACindex = asInfoKeys.Find("AC");
    ANindex = asInfoKeys.Find("AN");
    NSindex = asInfoKeys.Find("NS");

    int AN = 0; 
    int NS = 0;
    int ACs[3] = {0,0,0};
    for(int j=0; j < (int)vnSampleGenotypes.size(); ++j) {
      if ( vnSampleGenotypes[j] != 0xffff ) {
	AN += 2;
	++NS;
	++ACs[(vnSampleGenotypes[j] & 0x7f00) >> 8];
	if ( (vnSampleGenotypes[j] & 0x007f) != 0x007f ) {
	  ++ACs[(vnSampleGenotypes[j] & 0x007f)];
	}
      }
    }

    if ( ANindex >= 0 ) {
      asInfoValues[ANindex].printf("%d",AN);
    }
    else { //if ( AN > 0 ) {
      String tmp;
      tmp.printf("%d",AN);
      asInfoKeys.Add("AN");
      asInfoValues.Add(tmp);
    }

    if ( NSindex >= 0 ) {
      asInfoValues[NSindex].printf("%d",NS);
    }
    else { //if ( NS > 0 ) {
      String tmp;
      tmp.printf("%d",NS);
      asInfoKeys.Add("NS");
      asInfoValues.Add(tmp);
    }
    
    if ( ACindex >= 0 ) {
      if ( ACs[2] == 0 ) {
      asInfoValues[ACindex].printf("%d",ACs[1]);
    }
      else {
	asInfoValues[ACindex].printf("%d,%d",ACs[1],ACs[2]);
      }
    }
    else { // if ( ACs[0]+ACs[1]+ACs[2] > 0 ) {
      String tmp;
      if ( ACs[2] == 0 ) {
	tmp.printf("%d",ACs[1]);
      }
      else {
	tmp.printf("%d,%d",ACs[1],ACs[2]);
      }
      asInfoKeys.Add("AC");
      asInfoValues.Add(tmp);
    }
  }
}

void VcfMarker::setFormat(const String& s, bool upgrade) {
  // if upgrade is set, GT:GD:GQ are converted into GT:DP:GQ:PL
  if ( ( upgrade ) && ( s.Compare("GT:GD:GQ") == 0 ) ) {
    asFormatKeys.Clear();
    asFormatKeys.Add("GT");
    asFormatKeys.Add("DP");
    asFormatKeys.Add("GQ");
    asFormatKeys.Add("PL");
    GTindex = 0;
    GDindex = 1;
    GQindex = 2;
    DSindex = -1;
  }
  else {
    asFormatKeys.ReplaceColumns(s,':');
    GTindex = asFormatKeys.Find("GT");
    DSindex = asFormatKeys.Find("DS");
    GDindex = asFormatKeys.Find("DP");
    if ( GDindex < 0 ) 
      GDindex = asFormatKeys.Find("GD");
    GQindex = asFormatKeys.Find("GQ");
  }
}

void VcfMarker::setSampleSize(int newsize, bool parseGenotypes, bool parseDosages, bool parseValues) {
  //Logger::gLogger->error("VcfMarker::setSampleSize(%d)",newsize);
  nSampleSize = newsize;

  if ( parseValues ) {
    asSampleValues.Dimension(newsize * asFormatKeys.Length());
  }
  if ( parseGenotypes ) {
    vnSampleGenotypes.resize(newsize);
  }
  if ( parseDosages ) {
    vfSampleDosages.resize(newsize);
  }
}

void VcfMarker::setDosage(int sampleIndex, float dosage) {
  vfSampleDosages[sampleIndex] = dosage;
}

void VcfMarker::setGenotype(int sampleIndex, unsigned short genotype) {
  vnSampleGenotypes[sampleIndex] = genotype;
}

void VcfMarker::setSample(int sampleIndex, const String& sampleValue, bool parseGenotypes, bool parseDosages, bool parseValues, int minGD, int minGQ) {
  if ( !( parseValues || parseDosages || parseGenotypes ) ) {
    return;
  }
  //Logger::gLogger->error("VcfMarker::setSample(%d, %s, %d, %d, %d)",sampleIndex,sampleValue.c_str(),parseGenotypes,parseDosages,parseValues);

  if ( (sampleValue.Compare("./.") == 0) || ( sampleValue.Compare(".") == 0 ) ) {
    if ( parseValues ) {
      asSampleValues[asFormatKeys.Length() * sampleIndex] = "./.";
      for(int i=1; i < asFormatKeys.Length(); ++i) {
	asSampleValues[i + asFormatKeys.Length() * sampleIndex] = "";
      }
    }
    if ( parseGenotypes ) {
      vnSampleGenotypes[sampleIndex] = 0xffff;
    }
    if ( parseDosages ) {
      vfSampleDosages[sampleIndex] = -1; // missing
    }
  }
  else {
    tmpTokens.ReplaceColumns(sampleValue,':');
    if ( tmpTokens.Length() != asFormatKeys.Length() ) {
      throw VcfFileException("# values = %s do not match with # fields in FORMAT field = %d at sampleIndex = %d",sampleValue.c_str(),asFormatKeys.Length(),sampleIndex);
    }
    
    if ( parseValues ) {
      for(int i=0; i < tmpTokens.Length(); ++i) {
	asSampleValues[i + tmpTokens.Length() * sampleIndex] = tmpTokens[i];
      }
    }
    if ( parseGenotypes ) {
      if ( ( GTindex < 0 ) || ( tmpTokens[GTindex][0] == '.' ) || ( ( tmpTokens[GTindex].Length() > 2 ) && ( tmpTokens[GTindex][2] == '.' ) ) ) {
	// missing - 0xfff
	vnSampleGenotypes[sampleIndex] = 0xffff;
      }
      else {
	int GD = INT_MAX;
	if ( ( minGD > 0 ) && ( GDindex >= 0 ) ) {
	  GD = tmpTokens[GDindex].AsInteger();
	}

	int GQ = INT_MAX;
	if ( ( minGQ > 0 ) && ( GQindex >= 0 ) ) {
	  GQ = tmpTokens[GQindex].AsInteger();
	}

	//fprintf(stderr,"GD=%d,GQ=%d,GDindex=%d,GQindex=%d,minGD=%d,minGQ=%d\n",GD,GQ,GDindex,GQindex,minGD,minGQ);

	if ( ( GD < minGD ) || ( GQ < minGQ ) ) {
	  vnSampleGenotypes[sampleIndex] = 0xffff;
	  asSampleValues[asFormatKeys.Length() * sampleIndex] = "./.";
	  for(int i=1; i < asFormatKeys.Length(); ++i) {
	    asSampleValues[i + asFormatKeys.Length() * sampleIndex] = "";
	  }
	}
	else {
	  int sepPos = tmpTokens[GTindex].Find('|');
	  bool phased = false;
	  
	  if ( sepPos >= 0 ) {
	    phased = true;
	  }
	  else {
	    sepPos = tmpTokens[GTindex].Find('/');
	    phased = false;
	    if ( sepPos < 0 ) {
	      // allow haploid only for non-autosomal chromosomes
	      if ( VcfHelper::chromName2Num(sChrom) < 23 ) {
		throw VcfFileException("Cannot parse the genotype field " + tmpTokens[GTindex] + " at chromosome " + sChrom);
	      }
	    }
	  }
	  
	  if ( sepPos >= 0 ) {
	    tmpTokens[GTindex][sepPos] = '\0';
	  }
	  //Logger::gLogger->writeLog("GTindex = %d, %s, %s",GTindex, tmpTokens[GTindex].c_str(), tmpTokens[GTindex].c_str()+sepPos+1);

	  //Logger::gLogger->writeLog("sampleIndex = %d",sampleIndex);
	  
	  int n1 = atoi(tmpTokens[GTindex].c_str());
	  int n2 = (sepPos < 0) ? 0x00ff : atoi(tmpTokens[GTindex].c_str() + sepPos+1); // set second allele as missing if haploid
	  
	  if ( ( !phased ) && ( n1 > n2 ) ) {
	    int tmp = n1;
	    n1 = n2;
	    n2 = tmp;
	  }
	  
	  unsigned short g = (unsigned short)( ((n1 & 0x00ff) << 8) | (n2 & 0x00ff) | ((phased & 0x0001) << 15) );
	  
	  vnSampleGenotypes[sampleIndex] = g;
	}
      }
    }
    if ( parseDosages ) {
      if ( ( DSindex < 0 ) || ( tmpTokens[DSindex][0] == '.' ) ) {
	vfSampleDosages[sampleIndex] = -1; // missing
      }
      else {
	vfSampleDosages[sampleIndex] = static_cast<float>(atof(tmpTokens[DSindex].c_str()));
      }
    }
  }
}

bool BedFile::iterateMarker() {
  // read a marker information from BIM file
  if ( line.ReadLine(iBimFile) > 0 ) {
    ++nNumLines;
  }
  else {
    return false;
  }
  lineTokens.ReplaceTokens(line, " \t\r\n");

  VcfMarker* pMarker;

  if ( nBuffers == 0 ) {
    pMarker = new VcfMarker;
    vpVcfMarkers.push_back(pMarker);
    ++nNumMarkers;
    ++nHead;
  }
  else {
    // make a circular list with constant size nBuffer
    pMarker = vpVcfMarkers[nHead];
    nHead = (nHead+1) % nBuffers;
    ++nNumMarkers;
  }

  if ( lineTokens.Length() < 6 ) {
    throw VcfFileException("BIM file has only %d columns - %s",lineTokens.Length(),lineTokens[0].c_str());
  }
  
  try {

    String strChr = lineTokens[0];
    if ( strChr.Compare("23") == 0 ) {
      strChr = "X";
    }
    else if ( strChr.Compare("24") == 0 ) {
      strChr = "Y";
    }
    else if ( strChr.Compare("25") == 0 ) {
      strChr = "XY";
    }
    else if ( strChr.Compare("26") == 0 ) {
      strChr = "MT";
    }

    // try to resolve different conventions for chromosome names by using
    // (1) original chr
    // (2) strChr : 23 -> "X", 24 -> "Y", 25 -> "XY", 26 -> "MT"
    // (3) "chr"+ chr
    // (4) "chr"+ strChr
    
    pMarker->setPos(lineTokens[3]);

    genomeIndex_t markerIndex = genomeSequence.getGenomePosition( lineTokens[0].c_str(), pMarker->nPos );
    if ( markerIndex == INVALID_GENOME_INDEX ) {
      // routines specific to BED file format
      markerIndex = genomeSequence.getGenomePosition( strChr.c_str(), pMarker->nPos );

      if ( markerIndex == INVALID_GENOME_INDEX ) {
	strChr = String("chr") + strChr;
	markerIndex = genomeSequence.getGenomePosition( strChr.c_str(), pMarker->nPos );
	if ( markerIndex == INVALID_GENOME_INDEX ) {
	  throw VcfFileException("Cannot parse chromosome name "+lineTokens[0]+" in BIM file");
	}
      }
    }

    pMarker->setChrom(strChr);
    pMarker->setID(lineTokens[1]);

    // determine refBase/altBase
    // if refBase matches a1 or a2, use it
    // if none matches and either is zero
    // do not allow for triallelic sites
    char refBase = genomeSequence[markerIndex];
    char a1 = lineTokens[4][0];
    char a2 = lineTokens[5][0];
    char altBase = determineAltBase(refBase, a1,a2);

    pMarker->setRef(String(refBase));
    pMarker->setAlts(String(altBase));

    // read genotypes
    ifread(iFile, pBedBuffer, nBytes);

    pMarker->setSampleSize((int)vpVcfInds.size(),bParseGenotypes,bParseDosages,bParseValues);
    int AC = 0, AN = 0, NS = 0;

    for(int i=0; i < (int)vpVcfInds.size(); ++i) {
      char g = (pBedBuffer[i/4] >> ((i % 4)*2)) & 0x03;

      //Logger::gLogger->writeLog("i = %d, g = %d",i,(int)g);
      switch(g) {
      case 0: // 0/0 -> 0x00
	if ( bRefIsAllele1 ) {
	  pMarker->setGenotype(i,0x0000); // 0/0
	}
	else {
	  pMarker->setGenotype(i,0x0101); // 0/0
	  AC += 2;
	}
	AN += 2;
	++NS;
	break;
      case 2: // 0/1
	pMarker->setGenotype(i,0x0001); // 0/1
	++AC;
	AN += 2;
	++NS;
	break;
      case 1: // Missing
	pMarker->setGenotype(i,0xffff); // missing
	break;
      case 3:
	if ( bRefIsAllele1 ) {
	  pMarker->setGenotype(i,0x0101); // 1/1
	  AC += 2;
	}
	else {
	  pMarker->setGenotype(i,0x0000); // 0/0
	}
	AN += 2;
	++NS;
	break;
      default:
	throw VcfFileException("Error in parsing genotypes");
      }
    }

    pMarker->setQual("100");
    pMarker->setFilters("PASS");

    String info;
    if ( AN > 0 ) {
      info.printf("NS=%d;AC=%d;AN=%d;AF=%.6lf",NS,AC,AN,(double)AC/(double)AN);
    }
    else {
      info.printf("NS=%d;AC=%d;AN=%d",NS,AC,AN);
    }
    pMarker->setInfo(info, false, false);

    return true;
  }
  catch (VcfFileException exc) {
    // add the line number to the error message
    throw VcfFileException(exc.msg + " See line " + nNumLines + ".");
  }
}

char BedFile::determineAltBase(char refBase, char a1, char a2) {
    if ( a1 == '0' ) {
      if ( a2 == '0' ) {
	// altBase is unknown
	bRefIsAllele1 = true;
	return '.';
      }
      else if ( a2 == refBase ) {
	bRefIsAllele1 = true;
	return '.';
      }
      else { // no way to detect flip, and assume that strand is correct
	bRefIsAllele1 = false;
	return a2;
      }
    }
    else if ( a1 == refBase ) {
      bRefIsAllele1 = true;
      if ( a2 == '0' ) {
	return '.';
      }
      else {
	return a2;
      }
    }
    else if ( a2 == refBase ) {
      bRefIsAllele1 = false;
      if ( a1 == '0' ) {
	return '.';
      }
      else {
	return a1;
      }
    }
    else {
      if ( bAllowFlip ) {
	const char* fwd = "ACGT";
	const char* rev = "TGCA";
	char f1 = '.';
	char f2 = '.';

	for(int i=0; i < 4; ++i) {
	  if ( a1 == fwd[i] ) {
	    f1 = rev[i];
	  }
	  if ( a2 == fwd[i] ) {
	    f2 = rev[i];
	  }
	}

	if ( f1 == refBase ) {
	  bRefIsAllele1 = true;
	  return f2;
	}
	else if ( f2 == refBase ) {
	  bRefIsAllele1 = false;
	  return f1;
	}
      }
      //throw VcfFileException("The ref/alt allele does not match to reference genome");
      return '.';
    }
}

void VcfFile::printVCFHeader(IFILE oFile) {
  for(int i=0; i < getMetaCount(); ++i) {
    ifprintf(oFile,"##%s=%s\n",(const char*)getMetaKey(i), getMetaValue(i, "<na>").c_str());
  }
  ifprintf(oFile,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
  if ( ( getSampleCount() > 0 ) && ( !bSiteOnly ) ) {
    ifprintf(oFile,"\tFORMAT");
    for(int i=0; i < getSampleCount(); ++i) {
      ifprintf(oFile,"\t%s",vpVcfInds[i]->sIndID.c_str());
    }
  }
  ifprintf(oFile,"\n");
}

void VcfFile::printVCFHeaderSubset(IFILE oFile, std::vector<int>& subsetIndices) {
  //fprintf(stderr,"foo\n");
  for(int i=0; i < getMetaCount(); ++i) {
    ifprintf(oFile,"##%s=%s\n",getMetaKey(i).c_str(), getMetaValue(i, "<na>").c_str());
  }
  ifprintf(oFile,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
  ifprintf(oFile,"\tFORMAT");
  for(int j=0; j < (int)subsetIndices.size(); ++j) {
    int i = subsetIndices[j];
    ifprintf(oFile,"\t%s",vpVcfInds[i]->sIndID.c_str());
  }
  ifprintf(oFile,"\n");
}

void VcfFile::printBEDHeader(IFILE oBedFile, IFILE oFamFile) {
  for(int i=0; i < getSampleCount(); ++i) {
    if ( vpVcfInds[i]->sFamID.Length() == 0 ) {
      ifprintf(oFamFile,"%s",vpVcfInds[i]->sIndID.c_str());
    }
    else {
      ifprintf(oFamFile,"%s",vpVcfInds[i]->sFamID.c_str());
    }

    ifprintf(oFamFile,"\t%s",vpVcfInds[i]->sIndID.c_str());

    if ( vpVcfInds[i]->sFatID.Length() == 0 ) {
      ifprintf(oFamFile,"\t0");
    }
    else {
      ifprintf(oFamFile,"%s",vpVcfInds[i]->sFatID.c_str());
    }

    if ( vpVcfInds[i]->sMotID.Length() == 0 ) {
      ifprintf(oFamFile,"\t0");
    }
    else {
      ifprintf(oFamFile,"\t%s",vpVcfInds[i]->sMotID.c_str());
    }

    switch( vpVcfInds[i]->gender ) {
    case VcfInd::UNKNOWN:
      ifprintf(oFamFile,"\t0");
      break;
    case VcfInd::MALE:
      ifprintf(oFamFile,"\t1");
      break;
    case VcfInd::FEMALE:
      ifprintf(oFamFile,"\t2");
      break;
    default:
      throw VcfFileException("Unrecognized value for gender");
      break;
    }
    ifprintf(oFamFile,"\t-9\n");
  }
  ifclose(oFamFile);

  char magicNumbers[3] = {0x6c,0x1b,0x01};
  oBedFile->ifwrite(magicNumbers, 3);
}

void VcfHelper::printArrayJoin(IFILE oFile, const StringArray& arr, const char* sep, const char* empty, int start, int end) {
  for(int i=start; i < end; ++i) {
    if ( arr[i].Length() > 0 ) {
      if ( i > start ) {
	ifprintf(oFile,"%s",sep);
      }
      ifprintf(oFile,"%s",arr[i].c_str());
    }
  }
}

void VcfHelper::printArrayJoin(IFILE oFile, const StringArray& arr, const char* sep, const char* empty) {
  int len = arr.Length();
  if ( len == 0 ) {
    ifprintf(oFile,"%s",empty);
  }
  else if ( len == 1 ) {
    ifprintf(oFile,"%s",arr[0].c_str());
  }
  else {
    printArrayJoin(oFile,arr,sep,empty,0,len);
  }
}

void VcfHelper::printArrayDoubleJoin(IFILE oFile, const StringArray& arr1, const StringArray& arr2, const char* sep1, const char* sep2, const char* empty) {
  int len1 = arr1.Length();
  int len2 = arr2.Length();

  if ( len1 == len2 ) {
    printArrayDoubleJoin(oFile, arr1, arr2, sep1, sep2, empty, 0, len1);
  }
  else {
    throw VcfFileException("Inconsistency between arr1.Length() == %d and arr2.Length() == %d", len1, len2);
  }
}

void VcfHelper::printArrayDoubleJoin(IFILE oFile, const StringArray& arr1, const StringArray& arr2, const char* sep1, const char* sep2, const char* empty, int start, int end) {
  for(int i=start; i < end; ++i) {
    if ( i > start ) {
      ifprintf(oFile,"%s",sep1);
    }
    ifprintf(oFile,"%s%s%s",arr1[i].c_str(),sep2,arr2[i].c_str());
  }
}

// baseQ assume that 2 
void VcfMarker::printVCFMarker(IFILE oFile, bool siteOnly, int qGeno) {
  String line;

  ifprintf(oFile,"%s",sChrom.c_str());
  ifprintf(oFile,"\t%d",nPos);
  ifprintf(oFile,"\t%s",sID.c_str());
  ifprintf(oFile,"\t%s",sRef.c_str());

  if ( asAlts.Length() == 1 ) {
    ifprintf(oFile,"\t%s",asAlts[0].c_str());
  }
  else {
    ifprintf(oFile,"\t");
    VcfHelper::printArrayJoin(oFile, asAlts, ",", ".");
  }

  if ( fQual < 0 ) {
    ifprintf(oFile,"\t.");
  }
  else {
    ifprintf(oFile,"\t%.0f",fQual);
  }

  if ( asFilters.Length() == 1 ) {
    ifprintf(oFile,"\t%s",asFilters[0].c_str());
  }
  else {
    ifprintf(oFile,"\t");
    VcfHelper::printArrayJoin(oFile, asFilters, ";", "PASS");
  }

  ifprintf(oFile,"\t");
  if ( asInfoKeys.Length() == 0 ) {
    ifprintf(oFile,".");   
  }
  else {
    VcfHelper::printArrayDoubleJoin(oFile, asInfoKeys, asInfoValues, ";", "=", ".");
  }

  if ( !siteOnly ) {
    if ( asSampleValues.Length() > 0 ) {
      ifprintf(oFile,"\t");
      VcfHelper::printArrayJoin(oFile, asFormatKeys, ":", ".");
    
      for(int i=0; i < getSampleSize(); ++i) {
	ifprintf(oFile,"\t");
	VcfHelper::printArrayJoin(oFile, asSampleValues, ":", ".", i*asFormatKeys.Length(), (i+1)*asFormatKeys.Length());
      }
    }
    else if ( vnSampleGenotypes.size() > 0 ) { // converting from BED to VCF
      ifprintf(oFile,"\tGT");
      if ( qGeno > 0 ) ifprintf(oFile,":PL");
      for(int i=0; i < (int)vnSampleGenotypes.size(); ++i) {
	if ( vnSampleGenotypes[i] == 0xffff ) {
	  ifprintf(oFile,"\t./.");
	  if ( qGeno > 0 ) ifprintf(oFile,":0,0,0");
	}
	// special case for haploid
	else if ( (vnSampleGenotypes[i] & 0x00ff) == 0x00ff ) {
	  ifprintf(oFile,"\t%d",(vnSampleGenotypes[i] & 0x7f00) >> 8);
	  if ( qGeno > 0 ) ifprintf(oFile,":0,0,0");
	}
	else {
	  int g1 = (int)((vnSampleGenotypes[i] & 0x7f00) >> 8);
	  int g2 = (int)(vnSampleGenotypes[i] & 0x7f);
	  ifprintf(oFile,"\t%d/%d",g1,g2);
	  if ( qGeno > 0 ) {
	    ifprintf(oFile,":");
	    if ( g1 + g2 == 0 ) ifprintf(oFile,"%d,%d,%d",0,qGeno,qGeno+qGeno);
	    else if ( g1 + g2 == 1 ) ifprintf(oFile,"%d,%d,%d",qGeno,0,qGeno);
	    else ifprintf(oFile,"%d,%d,%d",qGeno+qGeno,qGeno,0);
	  }
	}
      }
    }
  }
  ifprintf(oFile,"\n");
}

void VcfMarker::printVCFMarkerSubset(IFILE oFile, std::vector<int>& subsetIndices, bool includeMono) {
  String line;
  //fprintf(stderr,"%s\n",(includeMono ? "foo" : "bar"));

  int ACindex = -1;
  int ANindex = -1;
  int NSindex = -1;

  ACindex = asInfoKeys.Find("AC");
  ANindex = asInfoKeys.Find("AN");
  NSindex = asInfoKeys.Find("NS");

  int AN = 0; 
  int NS = 0;
  int ACs[3] = {0,0,0};
  for(int j=0; j < (int)subsetIndices.size(); ++j) {
    int i = subsetIndices[j];
    if ( vnSampleGenotypes[i] != 0xffff ) { // not a missing genotype
      AN += 2;
      ++NS;
      ++ACs[(vnSampleGenotypes[i] & 0x7f00) >> 8];
      if ( (vnSampleGenotypes[i] & 0x007f) != 0x007f ) {
	++ACs[(vnSampleGenotypes[i] & 0x007f)];
      }
    }
  }

  if ( ( ACs[1] == 0 ) && ( ACs[2] == 0 ) && ( !includeMono ) ) {
    return;
  }

  if ( ANindex >= 0 ) {
    asInfoValues[ANindex].printf("%d",AN);
  }
  else { //if ( AN > 0 ) {
    String tmp;
    tmp.printf("%d",AN);
    asInfoKeys.Add("AN");
    asInfoValues.Add(tmp);
  }

  if ( NSindex >= 0 ) {
    asInfoValues[NSindex].printf("%d",NS);
  }
  else { //if ( NS > 0 ) {
    String tmp;
    tmp.printf("%d",NS);
    asInfoKeys.Add("NS");
    asInfoValues.Add(tmp);
  }
  
  if ( ACindex >= 0 ) {
    if ( ACs[2] == 0 ) {
      asInfoValues[ACindex].printf("%d",ACs[1]);
    }
    else {
      asInfoValues[ACindex].printf("%d,%d",ACs[1],ACs[2]);
    }
  }
  else { // if ( ACs[0]+ACs[1]+ACs[2] > 0 ) {
    String tmp;
    if ( ACs[2] == 0 ) {
      tmp.printf("%d",ACs[1]);
    }
    else {
      tmp.printf("%d,%d",ACs[1],ACs[2]);
    }
    asInfoKeys.Add("AC");
    asInfoValues.Add(tmp);
  }

  ifprintf(oFile,"%s",sChrom.c_str());
  ifprintf(oFile,"\t%d",nPos);
  ifprintf(oFile,"\t%s",sID.c_str());
  ifprintf(oFile,"\t%s",sRef.c_str());

  if ( asAlts.Length() == 1 ) {
    ifprintf(oFile,"\t%s",asAlts[0].c_str());
  }
  else {
    ifprintf(oFile,"\t");
    VcfHelper::printArrayJoin(oFile, asAlts, ",", ".");
  }

  if ( fQual < 0 ) {
    ifprintf(oFile,"\t.");
  }
  else {
    ifprintf(oFile,"\t%.0f",fQual);
  }

  if ( asFilters.Length() == 1 ) {
    ifprintf(oFile,"\t%s",asFilters[0].c_str());
  }
  else {
    ifprintf(oFile,"\t");
    VcfHelper::printArrayJoin(oFile, asFilters, ";", "PASS");
  }

  ifprintf(oFile,"\t");
  if ( asInfoKeys.Length() == 0 ) {
    ifprintf(oFile,".");   
  }
  else {
    VcfHelper::printArrayDoubleJoin(oFile, asInfoKeys, asInfoValues, ";", "=", ".");
  }

  if ( asSampleValues.Length() > 0 ) {
    ifprintf(oFile,"\t");
    VcfHelper::printArrayJoin(oFile, asFormatKeys, ":", ".");
    
    for(int j=0; j < (int)subsetIndices.size(); ++j) {
      int i = subsetIndices[j];
      ifprintf(oFile,"\t");
      VcfHelper::printArrayJoin(oFile, asSampleValues, ":", ".", i*asFormatKeys.Length(), (i+1)*asFormatKeys.Length());
    }
  }
  else if ( vnSampleGenotypes.size() > 0 ) {
    ifprintf(oFile,"\tGT",line.c_str());
    for(int j=0; j < (int)subsetIndices.size(); ++j) {
      int i = subsetIndices[j];
      if ( vnSampleGenotypes[i] == 0xffff ) {
	ifprintf(oFile,"\t./.");
      }
      else if ( (vnSampleGenotypes[i] & 0x00ff) == 0x00ff ) {
	ifprintf(oFile,"\t%d",(vnSampleGenotypes[i] & 0x7f00) >> 8);
      }
      else {
	ifprintf(oFile,"\t%d/%d",((vnSampleGenotypes[i] & 0x7f00) >> 8),(vnSampleGenotypes[i] & 0x007f));
      }
    }
  }
  ifprintf(oFile,"\n");
}

void VcfMarker::printBEDMarker(IFILE oBedFile, IFILE oBimFile, bool siteOnly) {
  if ( sID.Compare(".") == 0 ) {
    ifprintf(oBimFile,"%d\t%s:%d\t0\t%d\t%s\t%s\n",VcfHelper::chromName2Num(sChrom),sChrom.c_str(),nPos,nPos,sRef.c_str(),asAlts[0].c_str());
  }
  else {
    ifprintf(oBimFile,"%d\t%s\t0\t%d\t%s\t%s\n",VcfHelper::chromName2Num(sChrom),sID.c_str(),nPos,sRef.c_str(),asAlts[0].c_str());
  }

  int nBytes = (getSampleSize()+3)/4;
  char* genos = new char[nBytes]();
  int g1, g2;

  for(int i=0; i < getSampleSize(); ++i) {
    g1 = (vnSampleGenotypes[i] & 0x007f);
    g2 = ((vnSampleGenotypes[i] & 0x7f00) >> 8);

    // if haploid, set them as heterozygote
    // note that this procedure is irreversible for now
    // (not possile to convert BED to haploids)
    if ( g1 == 0x007f ) 
      g1 = 0;

    switch( g1 + g2 ) {
    case 0:
      genos[i/4] |= (0x0 << ((i%4)*2));
      break;
    case 1:
      genos[i/4] |= (0x2 << ((i%4)*2));
      break;
    case 2:
      genos[i/4] |= (0x3 << ((i%4)*2));
      break;
    default:
      genos[i/4] |= (0x1 << ((i%4)*2));
      break;
    }
  }
  oBedFile->ifwrite(genos,nBytes);
  delete [] genos;
}
