//
//  main.cpp
//  buildGlobalISA
//
//  Created by mark enstrom on 9/10/18.
//  Copyright © 2018 mark enstrom. All rights reserved.
//

#include <iostream>
#include "readAlignment.hpp"
#include "NWAlign.hpp"

int main(int argc, const char * argv[]) {
    if (argc != 2) {
        cout << "Usage: buildBlobalISA matafile\n";
        exit(0);
    }
    string metaFile = argv[1];
    //
    // read metafile from param 1
    //
    ifstream inFile;
    inFile.open(metaFile);
    if (!inFile) {
        cout << "Can't open metafile " << metaFile << endl;
        exit(0);
    }
    string subject = "";
    string projectPath = "";
    string s = "";
    //
    // read metafile
    //
    // subject:Z09132
    // path:<some path>
    // file1
    // file2 ...
    //
    vector<string> files = vector<string>();
    vector<ReadAlignment*> reads = vector<ReadAlignment*>();
    do {
        getline(inFile, s);
        if (!inFile) break;
        if (s[0] == '#') continue;
        size_t loc = s.find("subject:");
        if (loc != string::npos) {
            subject = s.substr(loc+8);
            continue;
        }
        loc = s.find("path:");
        if (loc != string::npos) {
            projectPath = s.substr(loc+5);
            continue;
        }
        files.push_back(s);
    } while (true);
    
    if (subject == "") {
        cout << "Error, could not read subject from metafile \n";
        exit(0);
    }
    if (projectPath == "") {
        cout << "Error, could not read path from metafile\n";
        exit(0);
    }
    //
    // read all alignment files
    //
    int index = 0;
    for (auto f:files) {
        ReadAlignment *pra = new ReadAlignment(index++,projectPath + "blat_files/" + f,f);
        bool bStatus = pra->read();
        if (!bStatus) {
            exit(0);
        }
        reads.push_back(pra);
    }
    CombinedInserts ci = CombinedInserts();
    //
    // combine alignments with matching refseq
    // then combine alignments with very similar refseq
    //
    bool bFirst = true;
    for (auto pra:reads) {
        if (bFirst) {
            ci.initAlignments(pra->_vAlign);
            bFirst = false;
        } else {
            cout << "add Alignments " << pra->_inputFile << "\n";
            ci.addAlignments(pra->_vAlign);
        }
    }
    //
    // very high value 7
    //
    ci.mergeClose(subject,7);
    //
    // normalize
    //
    ci.normalize((int)files.size());
    ci.buildGlobalList(reads,subject,projectPath);
    //
    // output
    //
    ci.printAlignments(reads,(int)files.size(),subject,projectPath);
    return 0;
}
