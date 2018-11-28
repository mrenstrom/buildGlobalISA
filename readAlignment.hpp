//
//  readAlignment.hpp
//  buildGlobalISA
//
//  Created by mark enstrom on 9/10/18.
//  Copyright © 2018 mark enstrom. All rights reserved.
//

#ifndef readAlignment_hpp
#define readAlignment_hpp

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <vector>
#include <map>
#include <unordered_map>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <list>
#include <stdlib.h>
#include <assert.h>
#include <ctime>
#include <unistd.h>
#include <algorithm>

using namespace std;



#define KMER 8
//-----------------------------------------------------------------------------------------------
//
// alClass = alignment type: S - single, M - Multiple, U - unaligned, X - invalid
//
//
//-----------------------------------------------------------------------------------------------
class Alignment {
public:
    bool    _valid;
    int     _fileNumber;
    int     _index;
    string  _alClass;
    string  _align;
    int     _oCount;
    int     _nCount;
    double  _frac;
    double  _log2;
    double  _fracLog2;
    double  _r;
    int     _span;
    string  _refseq;
    string  _chr;
    int     _startPos;
    string  _dir;
    string  _globalID = "X";
    int     _globalIndex = 0;
    string  _id;
    string  _mergeID;
    unordered_map<std::string,size_t> _seqMap;

    Alignment(int i,string refseq){
        _fileNumber = i;
        _valid = true;
        _refseq = refseq;
        //
        // build seq map (8)
        //
        size_t l = _refseq.length() - KMER;
        if (l > (40-KMER)) l = (40-KMER);
        for (int i = 0; i < (l-KMER);i++) {
            string s = _refseq.substr(i,KMER);
            ++_seqMap[s];
        }
    };
    
    bool similar(string s2) {
        for (int i = 0; i < (s2.length()-KMER); i += KMER) {
            if  (_seqMap.find(s2.substr(i,KMER)) != _seqMap.end()) {
                return true;
            }
        }
        return false;
    }
};




typedef Alignment* PAlignment;
typedef vector<PAlignment> VecAlign;
typedef VecAlign* PVecAlign;
typedef unordered_map<string,PAlignment> AlignMap;
typedef unordered_map<string,PVecAlign> RefMap;




//---------------------------------------------------------------------------
//
// sort mechanism
//
//
//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
//
//---------------------------------------------------------------------------
class AlignGroup {
public:
    PVecAlign _pVecAlign;
    bool      _valid;
    int       _count;
    string    _refseq;
    string    _align;
    string    _chr;
    string    _dir;
    int       _startPos;
    string    _alClass;
    int       _groupID;
    
    AlignGroup(int i) {
        _pVecAlign = new VecAlign();
        _valid = true;
        _alClass = "U";
        _align = "X";
        _groupID = i;
    }
    void addAlignment(PAlignment p);
    
    void verifyCount() {
        _count = 0;
        for (auto pv:*_pVecAlign) {
            _count += pv->_oCount;
        }
    }
    
    string refseq() {
        return _refseq;
    }
    
    string align() {
        return _align;
    }
    
    string alClass() {
        return _alClass;
    }
    
    bool similar(string s) {
        return _pVecAlign->at(0)->similar(s);
    }
    
    PAlignment top() {
        //
        // look for top non "U"
        //
        for (auto p:*_pVecAlign) {
            if (p->_alClass != "U") {
                return p;
            }
        }
        return(_pVecAlign->at(0));
    }
};

typedef AlignGroup *PAlignGroup;
typedef vector<PAlignGroup> VecAlignGroup;
typedef VecAlignGroup *PVecAlignGroup;
typedef unordered_map<string,PAlignGroup> MergeMap;


class MergeEntry {
public:
    string     _cause;
    string     _individual;
    PAlignment _p0;
    PAlignment _p1;
    MergeEntry(string cause, string s,PAlignment p0,PAlignment p1) {
        _p0 = p0;
        _p1 = p1;
        _cause = cause;
        _individual = s;
    }
};

class MergeLog{
public:
    vector<MergeEntry*> _entries;
    MergeLog() {
        _entries = vector<MergeEntry*>();
    }
    
    void add(string c, PAlignment p0,PAlignment p1) {
        MergeEntry *pme = new MergeEntry(c,"Single",p0,p1);
        _entries.push_back(pme);
    }
    
    void add(string c, PAlignGroup pm0,PAlignGroup pm1) {
        string s = to_string(pm0->_groupID) + " " + to_string(pm1->_groupID);
        MergeEntry *pme = new MergeEntry(c,s,pm0->top(),pm1->top());
        _entries.push_back(pme);
    }
};


class ReadAlignment {
public:
    bool read();
    string _inputFile = "";
    string _shortName;
    int _fileNumber;
    VecAlign _vAlign = VecAlign();
    
    ReadAlignment(int i,string fullName,string name) {
        _fileNumber = i;
        _inputFile = fullName;
        _shortName =name;
    }
};


class CombinedInserts {
public:
    RefMap _refMap = RefMap();
    MergeMap _mergeMap = MergeMap();
    MergeLog _mergeLog = MergeLog();
    VecAlignGroup _VecAlignGroup = VecAlignGroup();
    AlignMap _globalMap = AlignMap();
    VecAlign _globalVec = VecAlign();
    vector<vector<PAlignment>*> _outputFiles = vector<vector<PAlignment>*>();
    
    void initAlignments(VecAlign vAlign);
    void addAlignments(VecAlign vAlign);
    void mergeClose(string subject,int errorThreshold);
    void normalize(int s);
    void buildGlobalList(vector<ReadAlignment*>  reads,string subject,string path);
    void printAlignments(vector<ReadAlignment*> ,int s,string subject, string path);
};




#endif /* readAlignment_hpp */
