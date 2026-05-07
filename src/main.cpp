
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include "utils.hpp"
void pairend(std::string fastq_1, std::string fastq_2, std::string out_file,
	size_t compress_level = 6, size_t min_overlap = 10, size_t max_mismatches = 5,
	bool compress_output = false, bool enable_oes_test = true,
	double alpha = 1.0, double beta = -1.0,
	int n_threads = 4, int batch_size = 1000000);
void demultiplex(std::string fastq_file, std::string primer_table,
	std::string identified_output, std::string unidentified_output,
	int n_threads = 4, bool with_tag = true,
	bool compress_output = false, int compress_level = 6,
	int primer_mismatch = 3);
void dereplicate( std::string input_fasta,
	const std::string output_fasta,
	int n_threads,
	bool compress_output = false,
	int compress_level = 6,
	size_t batch_size = 100000);

void denoise(std::string input_fasta,
	std::string output_fasta,
	size_t sequence_min_L = 10,
	size_t sequence_max_L = 250,
	size_t min_count = 1,
	double threshold = 0.05,
	int n_threads = 32,
	size_t max_diff = 1,
	size_t max_indel = 1,
	bool keep_internal = true,
	bool keep_previous_stat = false,
	bool compress_output = false,
	bool enable_avx2 = true,
	bool debug = false,
	int compress_level = 6);

void count_records( std::string& fasfile);
int main() {
	std::string  workspace = "/media/data/luoyuan/Projects/METAeDNA/ideal_test/";
	std::string f1, f2;
	//f1 = "/media/data/luoyuan/metabarcoding/HTS/RawData/YN_2209_trnL/20231111-LGH-YZH-P4-trnL/231109-A00212A/trnL-LGH-YZH-P4-trnL-LFK12761_L1_1.fq.gz";
	//f2 = "/media/data/luoyuan/metabarcoding/HTS/RawData/YN_2209_trnL/20231111-LGH-YZH-P4-trnL/231109-A00212A/trnL-LGH-YZH-P4-trnL-LFK12761_L1_2.fq.gz";
	f1 = "~/cpp_debug/simulation/metabarcoding_data_simulation/Outputs/grinder_teleo1/grinder_teleo1_R1.fastq.gz";
	f2 = "~/cpp_debug/simulation/metabarcoding_data_simulation/Outputs/grinder_teleo1/grinder_teleo1_R2.fastq.gz";
	// PEAR test ------------------------------------------
	 
	 //std::string demultiplex_table = "/media/data/luoyuan/Projects/METAeDNA/simulated_demultiplex_table.csv";
	 /*
	 printf("Testing pairend...\n");
	 
	
	 bool is_fq1 = isFastqFile(f1); bool is_fq2 = isFastqFile(f2);
	 bool is_gz1 = isGzFile(f1), is_gz2 = isGzFile(f2);
	 FILE* fp1 = nullptr; gzFile gz1 = nullptr;
	 FILE* fp2 = nullptr; gzFile gz2 = nullptr;
	 FastqRecord rec1, rec2,merg;
	 if (is_gz1) gz1 = gzopen(expandTilde(f1).c_str(), "rb");
	 else        fp1 = fopen(expandTilde(f1).c_str(), "rb");
	 if (is_gz2) gz2 = gzopen(expandTilde(f2).c_str(), "rb");
	 else        fp2 = fopen(expandTilde(f2).c_str(), "rb");
	 size_t i = 0;
	 base_probability base_pb;
	 std::vector<FastqRecord> rec_v;
	 while (i < 100) {
		 FastqRecord tempr;
		 readRec(is_gz1, fp1, gz1, tempr, is_fq1);
		 rec_v.push_back(tempr);
		 i++;
	 }
	 if (is_gz1)  gzclose(gz1);
	 else        fclose(fp1);
	 if (is_gz1) gz1 = gzopen(expandTilde(f1).c_str(), "rb");
	 else        fp1 = fopen(expandTilde(f1).c_str(), "rb");
	 base_pb = calculate_base_probabilities(rec_v);

	 readRec(is_gz1, fp1,gz1,rec1, is_fq1); readRec(is_gz2, fp2, gz2, rec2, is_fq2);
	 size_t ovl,mismatch;
	 double overlap_score;
	 find_best_overlap(rec1, rec2, base_pb, 10, 3, mismatch, ovl, overlap_score, merg);
	 printf("Overlap length: %zu, Mismatch count: %zu, Overlap score: %.4f\n", ovl, mismatch, overlap_score);
	 std::cout << merg.seq;
	 */
	 // pair end ---------------------------
	pairend(f1, f2, "/media/data/luoyuan/Downloads/test.fastq",6, 10, 3, false, true, 1, -1,  32, 100000);
	/*
	demultiplex("/media/data/luoyuan/cpp_debug/simulation/metabarcoding_data_simulation/Outputs/grinder_teleo1/pairend.fastq",
		"/media/data/luoyuan/metabarcoding/HTS/trnL_ngsfilter.txt",
		"/media/data/luoyuan/cpp_debug/simulation/metabarcoding_data_simulation/Outputs/grinder_teleo1/identified.fas",
		"/media/data/luoyuan/cpp_debug/simulation/metabarcoding_data_simulation/Outputs/grinder_teleo1/unidentified.fas",32);

	
	dereplicate("/media/data/luoyuan/cpp_debug/simulation/metabarcoding_data_simulation/Outputs/grinder_teleo1/identified.fas",
		"/media/data/luoyuan/cpp_debug/simulation/metabarcoding_data_simulation/Outputs/grinder_teleo1/unique_trnl.fas", 10, false, 6);
	//*/
	 /*
	 denoise("/media/data/luoyuan/cpp_debug/simulation/metabarcoding_data_simulation/Outputs/grinder_teleo1/unique.fas",
		 "/media/data/luoyuan/cpp_debug/simulation/metabarcoding_data_simulation/Outputs/grinder_teleo1/denoise_trnL.fas");
		 */
	//denoise("/media/data/luoyuan/cpp_debug/simulation/metabarcoding_data_simulation/Outputs/grinder_teleo1/unique_big.fas",
	//	"/media/data/luoyuan/cpp_debug/simulation/metabarcoding_data_simulation/Outputs/grinder_teleo1/denoised_big.fas");
	//COUT << count_records(workspace + "pairend.fas") << std::endl;
	//COUT << count_records(workspace + "unique_singleton.fas") << std::endl;
}
