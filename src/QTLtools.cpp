/*Copyright (C) 2015 Olivier Delaneau, Halit Ongen, Emmanouil T. Dermitzakis
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.*/

#include "mode_cis/cis_data.h"
#include "mode_trans/trans_data.h"
#include "mode_match/match_data.h"
#include "mode_fenrich/fenrich_data.h"
#include "mode_correct/correct_data.h"
#include "mode_rtc/rtc_data.h"
#include "mode_pca/pca_data.h"
#include "mode_genrich/genrich_data.h"
#include "mode_extract/extract_data.h"
#include "mode_ase/ase_data.h"
#include "mode_quan/quan_data.h"
#include "mode_union/union_data.h"
#include "mode_bamstat/bamstat_data.h"

void printModes(){
    vrb.ctitle("Usage:");
    vrb.print("  QTLtools [mode] [options]");
    vrb.print("  eg: QTLtools cis --help");
    vrb.ctitle("Available modes:");
    vrb.print("  bamstat   Calculate basic QC metrics for BAM/SAM");
    vrb.print("  match     Match VCF genotypes to BAM/SAM file");
    vrb.print("  pca       Calculate principal components for a BED/VCF/BCF file");
    vrb.print("  correct   Covariate correction of a BED file");
    vrb.print("  cis       cis QTL analysis");
    vrb.print("  trans     trans QTL analysis");
    vrb.print("  fenrich   Functional enrichment for QTLs");
    vrb.print("  genrich   GWAS enrichment for QTLs");
    vrb.print("  rtc       Regulatory Trait Concordance analysis");
    vrb.print("  rtc-union Find the union of QTLs");
    vrb.print("  extract   Data extraction mode");
    vrb.print("  quan      Quantification mode");
    vrb.print("  ase       Measure allelic imbalance at every het genotype");
}

int main(int argc, char ** argv) {

	//1. Start timing
	timer running_timer;

	//2. Open LOG file if necessary
	for (int a = 1 ; a < argc - 1 ; a ++) {
		if ((strcmp(argv[a], "--log") == 0) && !vrb.open_log(string(argv[a+1]))) vrb.error("Impossible to open log file!");
		if (strcmp(argv[a], "--silent") == 0) vrb.set_silent();
	}

	//3. Print header on screen
	vrb.ctitle("QTLtools");
	vrb.bullet("Authors : Olivier DELANEAU / Halit ONGEN / Emmanouil DERMITZAKIS");
	vrb.bullet("Contact : olivier.delaneau@gmail.com");
	vrb.bullet("Webpage : XXX");
	vrb.bullet("Version : " + string(QTLTOOLS_VERSION));
	vrb.bullet("Date    : " + running_timer.date());

	//4. Switch mode
	vector < string > args;
    if (argc < 2){
        printModes();
        vrb.error("Not enough options, check command line!");
    }
	for (int a = 2 ; a < argc ; a ++) args.push_back(string(argv[a]));

	//5.1. CIS mode
	if (strcmp(argv[1], "cis") == 0) cis_main(args);

	//5.2. TRANS mode
	else if (strcmp(argv[1], "trans") == 0) trans_main(args);

	//5.3. MATCH mode
	else if (strcmp(argv[1], "match") == 0) match_main(args);

	//5.4. FENRICH mode
	else if (strcmp(argv[1], "fenrich") == 0) fenrich_main(args);

	//5.5. GENRICH mode
	else if (strcmp(argv[1], "genrich") == 0) genrich_main(args);

	//5.6. CORRECT mode
	else if (strcmp(argv[1], "correct") == 0) correct_main(args);

	//5.7. RTC mode
	else if (strcmp(argv[1], "rtc") == 0) rtc_main(args);
    
    //5.8. PCA mode
    else if (strcmp(argv[1], "pca") == 0) pca_main(args);

    //5.9. EXTRACT mode
    else if (strcmp(argv[1], "extract") == 0) extract_main(args);
    
    //5.10. RTC-UNION mode
    else if (strcmp(argv[1], "rtc-union") == 0) union_main(args);

    //5.11. QUANTIFICATION mode
    else if (strcmp(argv[1], "quan") == 0) quan_main(args);
    
    //5.12. ASE mode
    else if (strcmp(argv[1], "ase") == 0) ase_main(args);

    //5.13. BAMSTAT mode
    else if (strcmp(argv[1], "bamstat") == 0) bamstat_main(args);

	//5.14. UNRECOGNIZED mode
    else{
        printModes();
        vrb.error("Unrecognized QTLtools mode!");
    }

	//5. Terminate
	vrb.title("Running time: " + stb.str(running_timer.abs_time()) + " seconds");
	vrb.close_log();
}