#ifndef __CSG_VCF_FILE_H_ 
#define __CSG_VCF_FILE_H_

//////////////////////////////////////////////////////////////////////////////
// vcfcooker/VcfFile.h 
// (c) 2010 Hyun Min Kang, Matthew Flickenger, Matthew Snyder, Paul Anderson
//          Tom Blackwell, Mary Kate Trost, and Goncalo Abecasis
// 
// This file is distributed as part of the vcfCooker source code package   
// and may not be redistributed in any form, without prior written    
// permission from the author. Permission is granted for you to       
// modify this file for your own personal use, but modified versions  
// must retain this copyright notice and must not be distributed.     
// 
// Permission is granted for you to use this file to compile vcfCooker    
// 
// All computer programs have bugs. Use this file at your own risk.   
//
// Saturday November 10th, 2010 
//
// This header file provides interface to read/write VCF and BED files
// The detailed description of VCF file format can be found at
// http://
// The detailed description of BED file format can be found at
// http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml

#include <vector>
#include <string>
#include <exception>
#include <errno.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>

#include "GenomeSequence.h"
#include "InputFile.h"
#include "StringBasics.h"
#include "StringArray.h"
//#include "Logger.h"

namespace libVcf
{
///////////////////////////////////////////////////////////
// VcfFileException class 
// Creates a exception with printf style outputs
//////////////////////////////////////////////////////////
class VcfFileException : public std::exception {
 public: 
  VcfFileException(const char* format, ... ) {
    va_list args;
    va_start (args, format);
    msg.printf("VcfFileException: ");
    msg.vprintf(format, args);
    va_end(args);
  }

  ~VcfFileException() throw() {}
  const char* what() const throw() { return msg.c_str(); }
  String msg;
};

////////////////////////////////////////////////////////////////////////////////////////
// VcfInd class 
// class for accessing each individual's info
// VCF file, by default, contains only individual IDs
// BED file provides a .fam file providing family information
// (not supported yet) In VCF file, if a separate .fam file is provided, then
//                     the VcfInd class can contain family information
///////////////////////////////////////////////////////////////////////////////////////
class VcfInd {
 public:
  String sIndID;
  String sFamID;
  String sFatID;
  String sMotID;
  enum {UNKNOWN, MALE, FEMALE} gender;

 VcfInd(const String& indID) : sIndID(indID) { gender = UNKNOWN; };
  VcfInd(const String& indID, const String& famID, const String& fatID, const String& motID, const String& gender);
};

class VcfHelper {
 public:
  ////////////////////////////////////////////////////////////////////////////////////////
  // static member functions
  ////////////////////////////////////////////////////////////////////////////////////////
  static std::vector<double> vPhred2Err;
  static StringArray asChromNames;
  static std::vector<int> vnChromNums;
  static int char2twobit[255];

  static int chromName2Num(const String& chr);
  static int compareGenomicPos(const String& chr1, int pos1, const String& chr2, int pos2);
  static uint64_t str2TwoBits(const char* s, int len);

  static bool initPhred2Error(int maxPhred = 255);
  static bool initChromNamesNums();
  static bool initChar2TwoBit();

  static void printArrayJoin(IFILE oFile, const StringArray& arr, const char* sep, const char* empty);
  static void printArrayJoin(IFILE oFile, const StringArray& arr, const char* sep, const char* empty, int start, int end);
  static void printArrayDoubleJoin(IFILE oFile, const StringArray& arr1, const StringArray& arr2, const char* sep1, const char* sep2, const char* empty);
  static void printArrayDoubleJoin(IFILE oFile, const StringArray& arr1, const StringArray& arr2, const char* sep1, const char* sep2, const char* empty, int start, int end);
};


////////////////////////////////////////////////////////////////////////////////////////
// VcfMarker class 
// class for accessing each marker info, including genotypes, dosages, and other fields
////////////////////////////////////////////////////////////////////////////////////////
class VcfMarker {
 public:
  ////////////////////////////////////////////////////////////////////////////////////////
  // core member variables
  ////////////////////////////////////////////////////////////////////////////////////////
  String sChrom; // chromosome name
  int nPos;      // 1-based chromosomal coordinates
  String sID;    // Individual IDs of each variant
  String sRef;   // Strings representing reference Base
  StringArray asAlts; // Array of non-reference alleles
  float fQual;   // QUAL field (numeric, -1 if '.')
  StringArray asFilters;   // Array of filter elements
  StringArray asInfoKeys;   // Keys in the INFO field
  StringArray asInfoValues; // Values in the INFO field
  StringArray asFormatKeys; // Keys in the FORMAT fields
  StringArray asSampleValues; // Values (#FORMAT fields)*(#inds) of sample values
  std::vector<unsigned short> vnSampleGenotypes; // Genotypes by GT field
  std::vector<float> vfSampleDosages;   // Dosages by DS field

  ////////////////////////////////////////////////////////////////////////////////////////
  // core member functions (will stay as public)
  ////////////////////////////////////////////////////////////////////////////////////////
 VcfMarker() : GTindex(-1), DSindex(-1), GDindex(-1), GQindex(-1), nSampleSize(0), bPreserved(true) {}

  int getSampleSize() { return nSampleSize; }
  void setChrom(const String& s);
  void setPos(const String& s);
  void setID(const String& s);
  void setRef(const String& s);
  void setAlts(const String& s);
  void setQual(const String& s);
  void setFilters(const String& s);
  void setInfo(const String& s, bool upgrade, bool updateAC);
  void setFormat(const String& s, bool upgrade);
  void setSampleSize(int newsize, bool parseGenotypes, bool parseDosages, bool parseValue);
  void setDosage(int sampleIndex, float dosage);
  void setGenotype(int sampleIndex, unsigned short genotype);
  void setSample(int sampleIndex, const String& sampleValue, bool parseGenotypes, bool parseDosages, bool parseValues, int minGD, int minGQ);
  // print the marker info in VCF or BED format
  void printVCFMarker(IFILE oFile, bool siteOnly, int qGeno = 0);
  void printVCFMarkerSubset(IFILE oFile, std::vector<int>& subsetIndices, bool includeMono = false);
  void printBEDMarker(IFILE oBedFile, IFILE oBimFile, bool siteOnly);

  ////////////////////////////////////////////////////////////////////////////////////////
  // internal member variables
  ////////////////////////////////////////////////////////////////////////////////////////
  StringArray tmpTokens;  // temporary tokens for additional tokenization
  int GTindex;            // index of GT field
  int DSindex;            // index of DS field
  int GDindex;            // index of GD or DP field
  int GQindex;            // index of GQ field
  int nSampleSize;
  int bPreserved;         // indicate whether the INFO/FORMAT fields are preserved
};

class VcfFile {
 public:
  IFILE iFile;              // input/output file handle
  StringArray asMetaKeys;   // meta keys starting with '##' 
  StringArray asMetaValues; // values of meta-fields
  std::vector<VcfInd*> vpVcfInds; // individual info
  std::vector<VcfMarker*> vpVcfMarkers; // marker info (only buffered ones)
  int nBuffers;             // number of buffered lines
  int nNumMarkers;          // number of lines read so far
  bool bSiteOnly;       // read/write only site information (without genotypes)
  bool bParseGenotypes; // parse genotype values (GT)
  bool bParseDosages;   // parse dosage values (DS)
  bool bParseValues;    // parse all sample values as string
  bool bUpgrade;        // upgrade marker info
  bool bEOF;            // EOF marker flag
  int nMinGD;
  int nMinGQ;

  int nHead;                // internal variable to keep track of end of buffer
  String line;  // buffer line to store input line
  StringArray lineTokens; // line is tokenized to lineTokens

  VcfFile();
  virtual ~VcfFile();

  // opens a VCF file
  void openForRead(const char* filename, int nbuf = 1);
  // iterate a marker
  virtual bool iterateMarker();
  // get the last marker
  VcfMarker* getLastMarker() { return getLastMarker(0); }
  VcfMarker* getLastMarker(int nFromHead) { 
    if ( bEOF ) {
      return NULL;
    }
    //Logger::gLogger->writeLog("VcfFile::getLastMarker(%d), nHead=%d, nNumMarkers=%d, nBuffers=%d, accessing %d",nFromHead,nHead,nNumMarkers,nBuffers,(nHead-nFromHead-1) % nBuffers);
    else if ( nBuffers == 0 ) {
      return vpVcfMarkers[nHead-nFromHead-1];
    }
    else if ( nFromHead < nBuffers ) {
      /*
      if ( vpVcfMarkers[(nHead-nFromHead-1) % nBuffers] == NULL ) {
	throw VcfFileException("VcfFile::getLastMarker() - Null pointer exception, nHead=%d, nFromHead=%d, nBuffers=%d", nHead, nFromHead, nBuffers);
      }
      */
      return vpVcfMarkers[(nHead-nFromHead-1+nBuffers) % nBuffers]; 
    }
    else {
      throw VcfFileException("VcfFile::getLastMarker() - Index out of bound. nFromHead = %d, nBuffers = %d", nFromHead, nBuffers);
    }
  }

  int readLine();    // read a buffer of line
  int nNumLines;     // total number of lines

  // setting the characteristics of the reader : should be done before opening the file
  void reset();
  void setSiteOnly(bool siteOnly) { bSiteOnly = siteOnly; } // read only first 8 columns
  void setUpgrade(bool upgrade) { bUpgrade = upgrade; }     // convert from v3.3 (glfMultiples) to v4.0
  void setParseGenotypes(bool parseGenotypes) { bParseGenotypes = parseGenotypes; } // parse GT tag separately
  void setParseDosages(bool parseDosages) { bParseDosages = parseDosages; } // parse DS tag separately
  void setParseValues(bool parseValues) { bParseValues = parseValues; }     // parse individual's entry as strings. For example, if FORMAT is GT:DS:GL value is 0/1:1.000:30,0,32 then it is parsed as "0/1","1.000","30,0,32" .. 

  // parsing headers and meta lines
  void parseMeta();  // parse meta information
  void parseMetaLine();
  void upgradeMetaLines();
  void verifyMetaLines();   // check the sanity of meta lines
  void parseHeader();       
  void verifyHeaderLine(); 
  bool isMetaLine();
  bool isHeaderLine();
  int getMetaCount();
  int getMetaIndex(const String& meta);
  bool hasMetaKey(const String& meta);
  String getMetaKey(int offset);
  String getMetaValue(const String& meta);
  String getMetaValue(int offset);
  String getMetaValue(int offset, const String& ifmissing);
  String getMetaValue(const String& meta, const String& ifmissing);
  int getSampleCount();
  int getSampleIndex(const String& name);
  bool hasSample(String name);
  String getSampleID(int offset);

  void printVCFHeader(IFILE oFile);  // print headers in VCF format
  void printVCFHeaderSubset(IFILE oFile, std::vector<int>& subsetIndices);
  void printBEDHeader(IFILE oBedFile, IFILE oFamFile); // print headers in BED format
};

// BED file format 
class BedFile : public VcfFile {
 public:
  IFILE iBimFile;
  IFILE iFamFile;
  String sRefFile;
  bool bAllowFlip;
  bool bRefIsAllele1;
  char* pBedBuffer;
  int nBytes;
  GenomeSequence genomeSequence;

  void openForRead(const char* bfile, const char* reffile, int nbuf = 1);
  void openForRead(const char* bedFile, const char* bimFile, const char* famFile, const char* reffile, int nbuf = 1);

  BedFile();
  virtual ~BedFile();
  virtual bool iterateMarker();
  void setAllowFlip(bool b) { bAllowFlip = b; }
  char determineAltBase(char refBase, char a1, char a2);
};

}
#endif // __CSG_VCF_FILE_H_
