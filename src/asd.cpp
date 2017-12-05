/*
asd - a program to quickly calculate pairwise individual allele sharing distances
    Copyright (C) 2011  Zachary A Szpiech (szpiechz@umich.edu)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <pthread.h>
//#include <cstdio>
#include <cstdlib>
#include <cctype>
#include "gzstream.h"
#include "param_t.h"
#include "asd-cli.h"
#include "asd-data.h"
#include "errlog.h"
#include "asd-dist.h"
#include "pbar.h"

using namespace std;


int main(int argc, char *argv[])
{
#ifdef PTW32_STATIC_LIB
    pthread_win32_process_attach_np();
#endif

    param_t *params = getCLI(argc, argv);
    if (params == NULL) return 0;

    string outfile = params->getStringFlag(ARG_OUTFILE);
    LOG.init(outfile);
    LOG.log(getCommandLineString(argc, argv));
    LOG.log("Output file basename:", outfile);

    bool argerr = false;

    vector<string> comboFiles = params->getStringListFlag(ARG_COMBINE);
    if (comboFiles[0].compare(DEFAULT_COMBINE) != 0) {
        try {
            combine_partial_files(params);
        }
        catch (...) {
            return -1;
        }
        return 0;
    }

    string filename = params->getStringFlag(ARG_FILENAME);
    string tped_filename = params->getStringFlag(ARG_TPED_FILENAME);
    string tfam_filename = params->getStringFlag(ARG_TFAM_FILENAME);
    string vcf_filename = params->getStringFlag(ARG_VCF_FILENAME);
    bool STRU = (filename.compare(DEFAULT_FILENAME) != 0);
    bool TPED = (tped_filename.compare(DEFAULT_TPED_FILENAME) != 0);
    bool TFAM = (tfam_filename.compare(DEFAULT_TFAM_FILENAME) != 0);
    bool VCF = (vcf_filename.compare(DEFAULT_VCF_FILENAME) != 0);
    if (!check_file_type(STRU, TPED, TFAM, VCF)) {
        LOG.err("ERROR: Must specify either stru, vcf, or tped/tfam.");
        return -1;
    }
    if (STRU) LOG.log("Input stru file:", filename);
    if (TPED) LOG.log("Input tped file:", tped_filename);
    if (TFAM) LOG.log("Input tfam file:", tfam_filename);
    if (VCF) LOG.log("Input vcf file:", vcf_filename);
    bool GRM = params->getBoolFlag(ARG_GRM);
    bool BIALLELIC = params->getBoolFlag(ARG_BIALLELIC);
    if(GRM && !BIALLELIC){
        LOG.log("WARNING: --grm forces --biallelic.");
        BIALLELIC = true;
    }
    LOG.log("Calculate GRM:", GRM);
    LOG.log("Biallelic flag set (increases efficiency):", BIALLELIC);

    double MAF = params->getDoubleFlag(ARG_MAF);
    if(!check_maf(MAF)){
        LOG.err("ERROR: MAF must be >= 0 and <= 1.");
        return -1;
    }
    LOG.log("MAF filter:", MAF);

    int sort = params->getIntFlag(ARG_SORT);
    if (STRU) LOG.log("Individual ID column:", sort);

    int num_threads = params->getIntFlag(ARG_THREAD);
    if (!check_int_gt_0(num_threads)) {
        LOG.err("ERROR: Must have a positive number of threads.");
        return -1;
    }
    LOG.log("Number of threads:", num_threads);

    bool PRINT_PARTIAL = params->getBoolFlag(ARG_PARTIAL);
    bool PRINT_LOG = params->getBoolFlag(ARG_LOG);
    if (PRINT_PARTIAL && PRINT_LOG) {
        LOG.err("ERROR: Must choose only one of --partial, --log.");
        return -1;
    }
    if (!PRINT_PARTIAL) {
        LOG.log("Print allele sharing matrix:", !PRINT_PARTIAL);
        LOG.log("Log transformed:", PRINT_LOG);
    }
    else LOG.log("Print partial allele sharing matrix:", PRINT_PARTIAL);

    bool CALC_ALL_IBS = params->getBoolFlag(ARG_CALC_IBS);
    if(CALC_ALL_IBS && GRM){
        LOG.err("ERROR: --grm and --ibs incompatible.");
        return -1;
    }
    LOG.log("Output IBS matricies:", CALC_ALL_IBS);

    string STRU_MISSING = params->getStringFlag(ARG_STRU_MISSING);
    string TPED_MISSING = params->getStringFlag(ARG_TPED_MISSING);
    if (STRU) {
        LOG.log("STRU missing code:", STRU_MISSING);
    }
    else if (TPED) {
        LOG.log("TPED missing code:", TPED_MISSING);
    }

    bool PRINT_LONG = params->getBoolFlag(ARG_LONG_FORMAT);
    bool PRINT_LONG_IBS = params->getBoolFlag(ARG_IBS_LONG);
    if((PRINT_LONG || PRINT_LONG_IBS) && PRINT_PARTIAL){
        LOG.err("ERROR: --long and --long-ibs are not compatible with --partial.");
        return -1;
    }

    int nrows = 0;
    int ncols = 0;

    structure_data *data;
    if (STRU) {
        try {
            if (BIALLELIC) {
                data = readData_stru(filename, sort, nrows, ncols, MAF, STRU_MISSING);
            }
            else {
                data = readData_stru2(filename, sort, nrows, ncols, MAF, STRU_MISSING);
            }
        }
        catch (...) {
            return -1;
        }
    }
    else if (TPED) {
        try {
            if (BIALLELIC) {
                data = readData_tped_tfam(tped_filename, tfam_filename, nrows, ncols, MAF, TPED_MISSING);
            }
            else {
                data = readData_tped_tfam2(tped_filename, tfam_filename, nrows, ncols, MAF, TPED_MISSING);
            }
        }
        catch (...) {
            return -1;
        }
    }
    else if (VCF) {
        try {
            if (BIALLELIC) {
                data = readData_vcf(vcf_filename, nrows, ncols, MAF);
            }
            else {
                data = readData_vcf2(vcf_filename, nrows, ncols, MAF);
            }
        }
        catch (...) {
            return -1;
        }
    }

    LOG.log("Sample size:", nrows/2);
    LOG.log("Number of loci:", ncols);

    int nind = nrows / 2;

    Bar pbar;
    barInit(pbar, nind*num_threads, 78);

    init_storage(nind, CALC_ALL_IBS);

    unsigned int *NUM_PER_THREAD = make_thread_partition(num_threads, ncols);

    LOG.log("Starting calculations with", num_threads, false);
    LOG.log(" threads.");

    work_order_t *order;
    pthread_t *peer = new pthread_t[num_threads];
    unsigned int prev_index = 0;
    for (int i = 0; i < num_threads; i++)
    {
        order = new work_order_t;
        order->first_index = prev_index;
        order->last_index = prev_index + NUM_PER_THREAD[i];
        prev_index += NUM_PER_THREAD[i];
        order->stru_data = data;
        order->CALC_ALL_IBS = CALC_ALL_IBS;
        order->bar = &pbar;
        order->threads = num_threads;
        if(!GRM){
            if (BIALLELIC) {
                pthread_create(&(peer[i]),
                               NULL,
                               (void *(*)(void *))calc_pw_as_dist,
                               (void *)order);
            }
            else {
                pthread_create(&(peer[i]),
                               NULL,
                               (void *(*)(void *))calc_pw_as_dist2,
                               (void *)order);
            }
        }
        else{
            pthread_create(&(peer[i]),
                            NULL,
                            (void *(*)(void *))calc_grm,
                            (void *)order);
        }
    }

    for (int i = 0; i < num_threads; i++) pthread_join(peer[i], NULL);

    cerr << endl;

    finalize_calculations(nind, ncols, CALC_ALL_IBS, GRM);

    write_dist_matrix(outfile, nind, ncols, data->ind_names, PRINT_PARTIAL, PRINT_LOG, PRINT_LONG, GRM);

    if (CALC_ALL_IBS)
    {
        write_ibs_matrices(outfile, nind, ncols, data->ind_names, PRINT_PARTIAL, PRINT_LONG_IBS);
    }

    delete [] NUM_PER_THREAD;
    delete [] peer;

#ifdef PTW32_STATIC_LIB
    pthread_win32_process_detach_np();
#endif

    return 0;
}

