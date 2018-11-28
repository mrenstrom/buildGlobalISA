//
//  readAlignment.cpp
//  buildGlobalISA
//
//  Created by mark enstrom on 9/10/18.
//  Copyright © 2018 mark enstrom. All rights reserved.
//
#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>
#include <iomanip>
#include "readAlignment.hpp"
#include "NWAlign.hpp"
using namespace std;


struct myclassV {
    bool operator() (PAlignment pa1,PAlignment pa2) {
        return (pa1->_oCount > pa2->_oCount);
    }
} smallCompObj;

struct myclassM {
    bool operator() (PAlignGroup pm1,PAlignGroup pm2) {
        return (pm1->_count > pm2->_count);
    }
} smallCompObjM;
/*--------------------------------------------------------------------------------------------
 * add new alignment to group
 *
 *
 *--------------------------------------------------------------------------------------------*/
void AlignGroup::addAlignment(PAlignment p) {
    _pVecAlign->push_back(p);
    std::sort(_pVecAlign->begin(),_pVecAlign->end(),smallCompObj);
    int count = 0;
    //
    // find highest non-"U" alignment if there is one.
    //
    string al = "U";
    _refseq = _pVecAlign->at(0)->_refseq;
    // 
    // count-based totals
    //  
    int cS = 0;
    int cM = 0;
    int cU = 0;
    for (auto pv:*_pVecAlign) {
        //
        // full count
        //
        count += pv->_oCount;
        //
        // classify
        //
        if (pv->_alClass == "S") {
            cS+= pv->_oCount;
        } else if (pv->_alClass == "M") {
            cM+= pv->_oCount;
        } else if (pv->_alClass == "U") {
            cU+= pv->_oCount;
        }
    }
    if ((cS == 0) && (cM == 0)) {
        // unaligned
        _alClass = "U";
        _align   = "U";
        _chr = "U";
        _dir = "0";
        _startPos = 0;
    }
    //
    // can it be assigned 'S'?
    //
    else if (cS >= cM) {
        _alClass = "S";
        // find highest S
        for (auto pv:*_pVecAlign) {
            if (pv->_alClass == "S") {
                _align = pv->_align;
                _chr = pv->_chr;
                _dir = pv->_dir;
                _startPos = pv->_startPos;
                break;
            }
        }
    } else {
        _alClass = "M";
        // highest M
        for (auto pv:*_pVecAlign) {
            if (pv->_alClass == "M") {
                _align = pv->_align;
                _chr = pv->_chr;
                _dir = pv->_dir;
                _startPos = pv->_startPos;
                break;
            }
        }
    }
//    for (auto pv:*_pVecAlign) {
//        if ((al == "U") && (pv->_alClass != "U")) {
//            al = pv->_align;
//            _refseq = pv->_refseq;
//            _alClass = pv->_alClass;
//        }
//       count += pv->_oCount;
//    }
    _count = count;
}


/*--------------------------------------------------------------------------------------------
 *
 * read ISA alignment file...
 *
 *id                   count            alignStr    chromo     start     dir    span  rUnique       cUnique    refseq
 *MAPRIS:00003:571:97    571    M_chr1:31133282:-    M_chr1    31133282    -    97    0.975247524752    289    TCTCCTGCCTCAGCTTCCCAAGTACTGGGA
 *
 *--------------------------------------------------------------------------------------------*/
bool
ReadAlignment::read()
{
    ifstream inFile;
    inFile.open(_inputFile);
    if (!inFile) {
        cout << "Couldn't open fille " << _inputFile << endl;
        return false;
    }
    char line[1024];
    //
    // read first header line
    //
    inFile.getline(line, 1024);
    //
    // read rest of file
    //
    string s = "";
    string s_id;
    string s_class;
    string s_align;
    string s_index;
    string s_oCount;
    string s_nCount;
    string s_rUnique;
    string s_cUnique;
    string s_fracLog2;
    string s_span;
    string s_refseq;
    string s_chr;
    string s_startPos;
    string s_dir;
    string s_mergeID;
    
    do {
        getline(inFile, s_id, '\t');
        if (s_id.substr(0,1) == "#") {
            // comment, read rest of line
            getline(inFile, s_id);
            continue;
        }
        
        getline(inFile, s_oCount, '\t');
        getline(inFile, s_align, '\t');
        getline(inFile, s_chr,'\t');
        getline(inFile, s_startPos,'\t');


        getline(inFile, s_dir,'\t');
        getline(inFile, s_span,'\t');

        getline(inFile, s_rUnique, '\t');
        getline(inFile, s_cUnique, '\t');
        getline(inFile, s_refseq);
        if (!inFile) break;
        // 
        // check for refseq < 30...this would
        // only happen with a failure to find consensus sequence
        //
        if (s_refseq.length() < 30) {
            cout << "error in alignment : refseq length = " << s_refseq.length() << "\n";
            cout << _inputFile << endl;
            cout << s_id << " " << s_refseq << endl;
            //return false;
            continue;
        }
        
        Alignment *pal = new Alignment(_fileNumber,s_refseq);
        pal->_id       = s_id;
        pal->_mergeID  = s_mergeID;
        pal->_alClass  = "S";
        pal->_align    = s_align;
        pal->_oCount   = stoi(s_oCount);
        pal->_nCount   = 0;
        pal->_frac     = 0.0;
        pal->_log2     = 0.0;
        pal->_fracLog2 = 0.0;
        pal->_r        = stof(s_rUnique);
        if (pal->_r > 0.75) {
          pal->_alClass = "M";
        }
        pal->_span     = stoi(s_span);
        pal->_chr      = s_chr;
        pal->_startPos = stoi(s_startPos);
        pal->_dir      = s_dir;

        //cout << pal->_id << endl;

        //
        // make sure for now refseq = 30
        //
        if (pal->_refseq.length() > 30) {
            pal->_refseq = pal->_refseq.substr(0,30);
        }
        //
        // add single or multu-alignment or unaligned
        // just not "X" (already merged)
        //
        if ((pal->_alClass == "S") || (pal->_alClass == "M") || (pal->_alClass == "U")) {
            _vAlign.push_back(pal);
        }
        if (!inFile.good()) break;
    } while (true);
    
    return true;
}

/*--------------------------------------------------------------------------------------------
 *
 *
 *--------------------------------------------------------------------------------------------*/

bool alignmentMatch(AlignGroup *p1, AlignGroup* p2)
{
    if (p1->alClass() == "U") {
        return false;
    }
    if (p2->alClass() == "U") {
        return false;
    }
    //
    // is it multi?
    //
    if ((p1->alClass() == "M") || (p2->alClass() == "M")) {
        return true;
    }
    //
    // are alignments on same chr and dir, and close dist?
    //
    if ((p1->_chr) != (p2->_chr)) return false;
    if ((p1->_dir) != (p2->_dir)) return false;
    if (abs(p1->_startPos - p2->_startPos) > 50) return false;
    
    return true;
}

/*--------------------------------------------------------------------------------------------
 * initAlignments
 *  add first file to new alignment struncture
 *--------------------------------------------------------------------------------------------*/
void CombinedInserts::initAlignments(VecAlign vAlign) {
    cout << "Start initAlignmens\n";
    //
    // add to list
    //
    for (auto pal:vAlign) {
        try {
            // allow merges
            PVecAlign pval = _refMap.at(pal->_refseq);
            //cout << "-----Init Exact match--------------------------------------------\n";
            //cout << "    " << pal->_refseq << " " << pal->_nCount << endl;
            pval->push_back(pal);
            pal->_valid = false;
        } catch(std::out_of_range) {
            //
            // add new
            //
            PVecAlign pval = new VecAlign();
            pval->push_back(pal);
            _refMap[pal->_refseq] = pval;
            pal->_valid = false;
        }
    }
}
//---------------------------------------------------------------------------
// addAlignments
//   check for exact refseq match to existing alignment
//   check for very close refseq match to existing alignemnt
//
//---------------------------------------------------------------------------
void CombinedInserts::addAlignments(VecAlign vAlign)
{
    int exactMatch = 0;
    int newUnique = 0;
    int refMismatch = 0;
    int alignMatch = 0;
    //
    // add exact refseq matches to list
    //
    for (auto pal:vAlign) {
        //
        // if refseq already in map, add alignment to list of
        // alignments for exact match refseq
        //
        try {
            PVecAlign pval = _refMap.at(pal->_refseq);
            PAlignment p0 = pval->at(0);
            // do they match
            if ((pal->_chr != p0->_chr) ||
                (abs(pal->_startPos - p0->_startPos) > 10)) {
                refMismatch++;
//                cout << "exact refset alignment mismatch\n";
//                cout << p0->_id << "   " << p0->_alClass << "  " << p0->_refseq <<  endl;
//                cout << p0->_chr << endl;
//                cout << p0->_startPos << endl;
//
//
//                cout << pal->_id << "   " << pal->_alClass << "  " << pal->_refseq << endl;
//                cout << pal->_chr << endl;
//                cout << pal->_startPos << endl;
//
//                cout << "-------------------------------------\n";
            }
            
            
            pval->push_back(pal);
            pal->_valid = false;
            exactMatch++;
            _mergeLog.add("RefExact",p0,pal);
        } catch(std::out_of_range) {

        }
    }
    //
    // for sequences that did not have exact refseq match...(_valid = true)
    // look for exact alignment match
    //
    // !!! this is not sorted in order of size
    //
    for (auto pal2:vAlign) {
        if ((pal2->_valid) && (pal2->_alClass != "U")) {
            for (auto pr:_refMap) {
                PVecAlign pv = pr.second;
                Alignment* palRef = pv->at(0);
                string chrA = pal2->_chr;
                string chrB = palRef->_chr;
                
                if ((chrA == chrB) &&
                    (pal2->_dir  == palRef->_dir) &&
                    (abs(pal2->_startPos - palRef->_startPos) <= 50)) {
                    pal2->_valid = false;
                    pv->push_back(pal2);
                    _mergeLog.add("AlignMatch",palRef,pal2);
                    alignMatch++;
                    break;
                }
            }
        }
    }
    //
    // if sequence did not have refseq or alignment match with
    // sequences already in dictionary, add new
    //
    for (auto pal2:vAlign) {
        //
        // only add unaligned > 1
        //
        if ((pal2->_valid) && ((pal2->_alClass != "U") || ((pal2->_alClass == "U") && (pal2->_oCount > 1)))
            ) {
            try {
                //
                // already in dict...just added from the same file
                //
                PVecAlign pval = _refMap.at(pal2->_refseq);
                pval->push_back(pal2);
                pal2->_valid = false;
                _mergeLog.add("RefseqMatchSameFile",pval->at(0),pal2);
            } catch(std::out_of_range) {
                //
                // add new
                //
                PVecAlign pval = new VecAlign();
                pval->push_back(pal2);
                _refMap[pal2->_refseq] = pval;
                pal2->_valid = false;
                newUnique++;
            }
        }
    }
    cout << "\n";
    cout << "exact refseq match  = " << exactMatch << endl;
    cout << "refseq aln mismatch = " << refMismatch << endl;
    cout << "align match         = " << alignMatch << endl;
    cout << "new unique          = " << newUnique << endl;
}
/*--------------------------------------------------------------------------------------------
 * CombinedInserts::mergeClose
 *
 *  bassed on REFSEQ, merge alignments with similar refseq and merge-able alignments
 *      S - S    same : merged
 *      S - S    defferent : not merged
 *      S - M    merged
 *      M - M    merged
 *      U - S/M  merged
 *      U - U    merged
 *--------------------------------------------------------------------------------------------*/

void CombinedInserts::mergeClose(string subject,int errorThreshold) {
    cout << "Merge close alignments\n";
    int combineMatch = 0;
    //
    // convert map to vector
    //
    int id = 0;
    for (auto pr: _refMap) {
        PVecAlign pva = pr.second;
        PAlignGroup pmg = new AlignGroup(id++);
        for (auto pv:*pva) {
            pmg->addAlignment(pv);
            _mergeLog.add("assignGroup " + to_string(pmg->_groupID),pv,pv);
        }
       _VecAlignGroup.push_back(pmg);
    }
    
    _refMap.clear();
    
    std::sort(_VecAlignGroup.begin(), _VecAlignGroup.end(), smallCompObjM);
    
    cout << "mergesize = " << _VecAlignGroup.size() << endl;
//
// uncomment for debug
//
//    for (auto pm:_VecAlignGroup) {
//        if (pm->refseq() == "CAGTATATGATTCAAAAAAAAAAAGCTATG") {
//            cout << "found\n";
//            for (auto p:*pm->_pVecAlign) {
//                cout << "    " << p->_refseq << " " << p->_align << " " << p->_oCount << endl;
//            }
//        }
//        if (pm->_count > 100) {
//            cout << pm->refseq() << " " << pm->align() << " " << pm->_count << " " << i++ << endl;
//        }
//    }
    
    //
    // from largest to smallest
    //
    for (int i = 0; i < _VecAlignGroup.size(); i++) {
        PAlignGroup pmA = _VecAlignGroup[i];
        
        if ((pmA->_valid) && (pmA->_count > 1)) {
            //cout << "Start Search " << i << " " << pmA->refseq() << " " << pmA->_count << endl;
            //
            // from smallest to large
            //
            for (int j = (int)_VecAlignGroup.size()-1; j > i; j--) {
                PAlignGroup pmB = _VecAlignGroup[j];
                
                if (pmB->_valid) {
                    bool show = false;
                   
                    if (
                         (pmB->refseq().substr(0,9)== "GTGTGGTGG") && 
                         (pmA->refseq().substr(0,9)== "GTGTGGTGG")){ 
                        cout << pmA->refseq() << endl;
                        cout << pmB->refseq() << endl;
                        NWAlign nw = NWAlign();
                        int dist = nw.align(pmA->refseq(),pmB->refseq());
                        cout << "i = " << i << " j = " << j << " x dist = " << dist << endl;
                        show = true;
                    }
                    
                    if (pmA->similar(pmB->refseq())) {
                        NWAlign nw = NWAlign();
                        int dist = nw.align(pmA->refseq(),pmB->refseq());
                        if ((dist <= errorThreshold)) {
                            if (show) {
                               cout << "Dist ok " << dist << endl;
                            }
                            if (alignmentMatch(pmA,pmB)) {
                                if (show) cout << "Alignment match\n";
                                //
                                // add all alignments from B to A
                                //
                                for (auto pmB_element:*pmB->_pVecAlign) {
                                    pmA->addAlignment(pmB_element);
                                }
                                //
                                // mark B as invalid (merged)
                                //
                                pmB->_valid = false;
                                combineMatch++;
                                //
                                // log it
                                //
                                _mergeLog.add("RefseqClose",pmA, pmB);
                            } else {
                               if (show) cout << "Alignment MISMATCH\n";
                            }
                        }
                    }
                }
            }
        }
    }
    //
    // get rid of merged alignments (marked invalid)
    // also get rid of unaligned sequences that couldn't be merged
    //
    VecAlignGroup tVec = VecAlignGroup();
    for (auto pm: _VecAlignGroup) {
        cout << pm->_align << " " << pm->_valid << " " << pm->alClass() << endl;
        if ((pm->_valid) && (pm->alClass() != "U")) {
            tVec.push_back(pm);
            cout << "add \n";
        } else {
            if (pm->alClass() == "U") {

               cout << "drop \n";
                _mergeLog.add("DropUnaligned " + to_string(pm->_groupID),pm->top(), pm->top());
            }
        }
    }
    _VecAlignGroup.clear();
    for (auto pm:tVec) {
        _VecAlignGroup.push_back(pm);
    }
    tVec.clear();
    //
    // verify count
    //
    for (auto pm:_VecAlignGroup) {
        pm->verifyCount();
    }
    //
    // re-sort
    //
    std::sort(_VecAlignGroup.begin(), _VecAlignGroup.end(), smallCompObjM);
    //
    // assign ID
    //
    int globalID = 1;
    for (auto pm:_VecAlignGroup) {
        char sg[1000];
        sprintf(sg,"%06i",globalID);
        string sGlobalID = subject + "_" + sg;
        for (auto pa: *pm->_pVecAlign) {
            pa->_globalID = sGlobalID;
            pa->_globalIndex = globalID;
        }
        _mergeMap[sGlobalID] = pm;
        globalID++;
    }
}

//---------------------------------------------------------------------------
// CombinedInserts::normalize
//
//
//
//---------------------------------------------------------------------------

void CombinedInserts::normalize(int s)
{
    //
    // output files...each file is a PVecAlign
    // list of files is array of pVecAlign
    //
    cout << "Build output file list \n";
    for (int i = 0; i < s; i++) {
        vector<PAlignment> *pvFile = new vector<PAlignment>();
        _outputFiles.push_back(pvFile);
    }
    //
    // for each entry in alignment map, go through all alignments and
    // assign to proper output file
    //

    for (auto pm:_VecAlignGroup) {
        PVecAlign pva = pm->_pVecAlign;
        for (auto pa:*pva) {
            vector<PAlignment>* pvFile = _outputFiles.at(pa->_fileNumber);
            pvFile->push_back(pa);
            //cout << pa->_globalID << " " << pa->_refseq << " " << pa->_align << " " << pa->_oCount << " " <<  pa->_fileNumber << endl;
        }
    }
    //
    // for each file....
    //      merge alignments with same ID
    //      calculate total count
    //      calculate normalized count (as if each file sequenced to depth 100,000)
    //
    for (auto pva:_outputFiles) {
        //
        // pOutVec is the list of alignments assigned to a single output file.
        // multiple files may have been combined into one file (sumAlignments)
        // and the same global alignments may map to a different genomic location
        //
        // combine same global entries with different genomic alignments
        //
        unordered_map<string, PAlignment> fileMap = unordered_map<string, PAlignment>();
        for (auto pa:*pva) {
            string key = pa->_globalID;
            try {
                PAlignment paOld = fileMap.at(key);
                // already exists, merge
                //cout << "-----Merge final Alignment--------------------------------------------\n";
                //cout << "    " << paOld->_globalID << " " << paOld->_align << endl;
                //cout << "    " << pa->_globalID << " " << pa->_align << endl;
                paOld->_oCount += pa->_oCount;
                //pa->_oCount = 0;
                //pa->_globalID  = "";
            } catch(std::out_of_range) {
                //
                // add new
                //
                fileMap[key] = pa;
            }
        }
        //
        // convert map back to vector
        //
        pva->clear();
        for (auto pr:fileMap) {
            pva->push_back(pr.second);
        }
        //
        // calc total count for each
        //
        int totalCount = 0;
        for (auto pa:*pva) {
            totalCount += pa->_oCount;
        }
        //
        // normalize
        //
        double normalize = 100000.0/double(totalCount);
        for (auto pa:*pva) {
            pa->_nCount = int(pa->_oCount * normalize);
            if (pa->_nCount < 1) pa->_nCount = 1;
            pa->_frac   = (double)pa->_oCount / totalCount;
        }
    }
}


//---------------------------------------------------------------------------
//
// build global list
//
//
//
//---------------------------------------------------------------------------
void CombinedInserts::buildGlobalList(vector<ReadAlignment*>  reads,string subject,string path)
{   //
    //
    //
    string s1 = path + "final_files/globalID.txt";
    ofstream myStats;
    myStats.open(s1);
    for (auto pm: _VecAlignGroup) {
        myStats << pm->top()->_globalID << "\t" << pm->_count << "\t" << pm->refseq() <<  "\t" << pm->alClass() << "\t" << pm->align() << endl;
    }
    myStats.close();
    //
    // full global alignments
    //
    s1 = path + "final_files/globalAlignments.txt";
    myStats.open(s1);
    myStats << "ID\tRefSeq\toCount\tnCount\talign\tsource\n";
    for (auto pm: _VecAlignGroup) {
        for (auto pa: *pm->_pVecAlign) {
            myStats << pa->_globalID << "\t";
            myStats << pa->_refseq << "\t";
            myStats << pa->_oCount << "\t";
            myStats << pa->_nCount << "\t";
            myStats << pa->_alClass << "\t";
            myStats << pa->_align << "\t";
            ReadAlignment * pra = reads.at(pa->_fileNumber);
            myStats << pra->_shortName << endl;
        }
    }
    myStats.close();
    //
    // merge list
    //
    s1 = path + "final_files/mergeLog.txt";
    myStats.open(s1);
    for (auto pme:_mergeLog._entries) {
        myStats << pme->_cause << "\t";
        myStats << pme->_individual << "\t";
        myStats <<  pme->_p0->_globalID << "\t";
        myStats << pme->_p0->_id << "\t";
        myStats << pme->_p1->_id << "\t";
        myStats << pme->_p1->_refseq << endl;
    }
    myStats.close();
}
//---------------------------------------------------------------------------
//
// build global list
//
// build each file list
//
//---------------------------------------------------------------------------

void CombinedInserts::printAlignments(vector<ReadAlignment*>  reads, int s,string subject, string path)
{
    //
    // write out files
    //
    cout << "write output files\n";
    int i = 0;
    for (auto pOutVec:_outputFiles) {
        //
        // build output file name
        //
        ofstream myStats;
        string s1 = path;
        ReadAlignment * pa = reads.at(i++);
        size_t last = pa->_inputFile.find_last_of('/');
        string fname = pa->_inputFile.substr(last+1);
        // remove .tran
        fname = fname.substr(0,fname.length()-5);
        fname = fname + ".tsv";
        
        s1 = s1 + "final_files/" + fname;
        myStats.open(s1);
        //
        // create header
        //
        myStats << "id"          << "\t";
        myStats << "GlobalIndex" << "\t";
        myStats << "oCount"      << "\t";
        myStats << "nCount"      << "\t";
        myStats << "frac"        << "\t";
        myStats << "log2"        << "\t";
        myStats << "fracLog2"    << "\t";
        myStats << "span"        << "\t";
        myStats << "refseq"      << "\t";
        myStats << "align"       << "\n";
        //
        // sort by nCount
        //
        std::sort(pOutVec->begin(),pOutVec->end(),smallCompObj);
        //
        // output each line
        //
        for (auto pa:*pOutVec) {
            if ((pa->_alClass == "S") || (pa->_alClass == "M") || (pa->_alClass == "U")) {
                myStats << pa->_globalID << "\t";
                myStats << pa->_globalIndex << "\t";
                myStats << pa->_oCount   << "\t";
                myStats << pa->_nCount   << "\t";
                myStats << pa->_frac     << "\t";
                myStats << pa->_log2     << "\t";
                myStats << pa->_fracLog2 << "\t";
                myStats << pa->_r        << "\t";
                myStats << pa->_refseq   << "\t";
                if (pa->_alClass == "S") {
                    myStats << pa->_align;
                } else if (pa->_alClass == "M") {
                    myStats << "M_" << pa->_align;
                } else if (pa->_alClass == "U") {
                    //
                    // get top align from group
                    //
                    PAlignGroup pm = _mergeMap[pa->_globalID];
                    myStats << pm->align();
                } else {
                    cout << "Error bad class\n";
                }
                myStats << endl;
            }
        }
        myStats.close();
    }
}





////---------------------------------------------------------------------------
////
//// build global list
////
////
////
////---------------------------------------------------------------------------
//void CombinedInserts::buildGlobalList(vector<ReadAlignment*>  reads,string subject,string path)
//{
//    //
//    // build global list from file lists
//    //
//    for (auto pva:_outputFiles) {
//        for (auto pa:*pva) {
//            string key = pa->_globalID;
//
//            try {
//                PAlignment paGlobal = _globalMap.at(key);
//                paGlobal->_nCount += pa->_nCount;
//
//            } catch(std::out_of_range) {
//                //
//                // add new
//                //
//                PAlignment paNew = new Alignment(0,pa->_refseq);
//                paNew->_globalID = pa->_globalID;
//                int nCount = pa->_nCount;
//                // individual file can add max of 5% of 100,000 to total count
//                if (nCount > 5000) nCount = 5000;
//                paNew->_nCount = nCount;
//                _globalMap[key] = paNew;
//            }
//        }
//    }
//    //
//    // convert global map to vector
//    //
//    for (auto pr:_globalMap) {
//        _globalVec.push_back(pr.second);
//    }
//    //
//    // sort global on nCount
//    //
//    std::sort(_globalVec.begin(),_globalVec.end(),smallCompObj);
//    //
//    // translate all files to new sorted global IDs
//    //
//    int globalID = 1;
//    string s1 = path + "final_files/globalID.txt";
//    ofstream myStats;
//    myStats.open(s1);
//    for (auto paGlobal: _globalVec) {
//        char sg[1000];
//        sprintf(sg,"%06i",globalID);
//        string sGlobalID = subject + "_" + sg;
//        for (auto pOutVec:_outputFiles) {
//            for (auto pa:*pOutVec) {
//                if (paGlobal->_globalID == pa->_globalID) {
//                    pa->_globalID = sGlobalID;
//                    pa->_globalIndex = globalID;
//                }
//            }
//        }
//        myStats << sGlobalID << "\t" << paGlobal->_nCount << "\t" << paGlobal->_refseq <<  endl;
//        globalID++;
//    }
//    myStats.close();
//    //
//    // full global alignments
//    //
//    s1 = path + "final_files/globalAlignments.txt";
//    myStats.open(s1);
//    myStats << "ID\tRefSeq\toCount\tnCount\talign\tsource\n";
//    for (auto paGlobal: _globalVec) {
//        string key = paGlobal->_refseq;
//
//        PVecAlign pva = _refMap.at(key);
//        for (auto pa:*pva) {
//            // skip if alignment is merged into parent due to duplicates in one file
//            if (pa->_globalID == "") continue;
//            myStats << setw(15) << pa->_globalID << "\t";
//            myStats << key << "\t";
//            myStats << setw(6) << pa->_oCount << "\t";
//            myStats << setw(6) << pa->_nCount << "\t";
//            myStats << pa->_align << "\t";
//            ReadAlignment * pra = reads.at(pa->_fileNumber);
//            myStats << pra->_shortName << endl;
//        }
//    }
//    myStats.close();
//}


