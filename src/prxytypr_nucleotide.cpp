#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "prxytypr_nucleotide.h"
#include <string.h>
#include "prxytypr_ansi_string.h"

char codons[64][4] = { "AAA",
"AAC",
"AAG",
"AAT",
"ACA",
"ACC",
"ACG",
"ACT",
"AGA",
"AGC",
"AGG",
"AGT",
"ATA",
"ATC",
"ATG",
"ATT",
"CAA",
"CAC",
"CAG",
"CAT",
"CCA",
"CCC",
"CCG",
"CCT",
"CGA",
"CGC",
"CGG",
"CGT",
"CTA",
"CTC",
"CTG",
"CTT",
"GAA",
"GAC",
"GAG",
"GAT",
"GCA",
"GCC",
"GCG",
"GCT",
"GGA",
"GGC",
"GGG",
"GGT",
"GTA",
"GTC",
"GTG",
"GTT",
"TAA",
"TAC",
"TAG",
"TAT",
"TCA",
"TCC",
"TCG",
"TCT",
"TGA",
"TGC",
"TGG",
"TGT",
"TTA",
"TTC",
"TTG",
"TTT" };

char AAs[64][5] = { "Lys",
"Asn",
"Lys",
"Asn",
"Thr",
"Thr",
"Thr",
"Thr",
"Arg",
"Ser",
"Arg",
"Ser",
"Ile",
"Ile",
"Met",
"Ile",
"Gln",
"His",
"Gln",
"His",
"Pro",
"Pro",
"Pro",
"Pro",
"Arg",
"Arg",
"Arg",
"Arg",
"Leu",
"Leu",
"Leu",
"Leu",
"Glu",
"Asp",
"Glu",
"Asp",
"Ala",
"Ala",
"Ala",
"Ala",
"Gly",
"Gly",
"Gly",
"Gly",
"Val",
"Val",
"Val",
"Val",
"STOP",
"Tyr",
"STOP",
"Tyr",
"Ser",
"Ser",
"Ser",
"Ser",
"STOP",
"Cys",
"Trp",
"Cys",
"Leu",
"Phe",
"Leu",
"Phe" };

//TTT	Phe		TCT	Ser		TAT	Tyr		TGT	Cys
//TTC	Phe		TCC	Ser		TAC	Tyr		TGC	Cys
//TTA	Leu		TCA	Ser		TAA	STOP	TGA	STOP
//TTG	Leu		TCG	Ser		TAG	STOP	TGG	Trp
//CTT	Leu		CCT	Pro		CAT	His		CGT	Arg
//CTC	Leu		CCC	Pro		CAC	His		CGC	Arg
//CTA	Leu		CCA	Pro		CAA	Gln		CGA	Arg
//CTG	Leu		CCG	Pro		CAG	Gln		CGG	Arg
//ATT	Ile		ACT	Thr		AAT	Asn		AGT	Ser
//ATC	Ile		ACC	Thr		AAC	Asn		AGC	Ser
//ATA	Ile		ACA	Thr		AAA	Lys		AGA	Arg
//ATG	Met		ACG	Thr		AAG	Lys		AGG	Arg
//GTT	Val		GCT	Ala		GAT	Asp		GGT	Gly
//GTC	Val		GCC	Ala		GAC	Asp		GGC	Gly
//GTA	Val		GCA	Ala		GAA	Glu		GGA	Gly
//GTG	Val		GCG	Ala		GAG	Glu		GGG	Gly 
char** copy_codons()
{
	 char codons[64][4] = { "AAA",
	"AAC",
	"AAG",
	"AAT",
	"ACA",
	"ACC",
	"ACG",
	"ACT",
	"AGA",
	"AGC",
	"AGG",
	"AGT",
	"ATA",
	"ATC",
	"ATG",
	"ATT",
	"CAA",
	"CAC",
	"CAG",
	"CAT",
	"CCA",
	"CCC",
	"CCG",
	"CCT",
	"CGA",
	"CGC",
	"CGG",
	"CGT",
	"CTA",
	"CTC",
	"CTG",
	"CTT",
	"GAA",
	"GAC",
	"GAG",
	"GAT",
	"GCA",
	"GCC",
	"GCG",
	"GCT",
	"GGA",
	"GGC",
	"GGG",
	"GGT",
	"GTA",
	"GTC",
	"GTG",
	"GTT",
	"TAA",
	"TAC",
	"TAG",
	"TAT",
	"TCA",
	"TCC",
	"TCG",
	"TCT",
	"TGA",
	"TGC",
	"TGG",
	"TGT",
	"TTA",
	"TTC",
	"TTG",
	"TTT" };

	char** codons_copy = new char*[64];
	for (int i_trip = 0; i_trip < 64; i_trip++)
	{
		codons_copy[i_trip] = t_string::copy_me_str(codons[i_trip]);
	} // AA_copy loop.

	return(codons_copy);
} // copy_codons loop.

char** copy_AAs()
{
	char AAs[64][5] = { "Lys",
	"Asn",
	"Lys",
	"Asn",
	"Thr",
	"Thr",
	"Thr",
	"Thr",
	"Arg",
	"Ser",
	"Arg",
	"Ser",
	"Ile",
	"Ile",
	"Met",
	"Ile",
	"Gln",
	"His",
	"Gln",
	"His",
	"Pro",
	"Pro",
	"Pro",
	"Pro",
	"Arg",
	"Arg",
	"Arg",
	"Arg",
	"Leu",
	"Leu",
	"Leu",
	"Leu",
	"Glu",
	"Asp",
	"Glu",
	"Asp",
	"Ala",
	"Ala",
	"Ala",
	"Ala",
	"Gly",
	"Gly",
	"Gly",
	"Gly",
	"Val",
	"Val",
	"Val",
	"Val",
	"STOP",
	"Tyr",
	"STOP",
	"Tyr",
	"Ser",
	"Ser",
	"Ser",
	"Ser",
	"STOP",
	"Cys",
	"Trp",
	"Cys",
	"Leu",
	"Phe",
	"Leu",
	"Phe" };

	char** AA_copy = new char*[64];
	for (int i_trip = 0; i_trip < 64; i_trip++)
	{
		AA_copy[i_trip] = t_string::copy_me_str(AAs[i_trip]);
	} // AA_copy loop.

	return(AA_copy);
} // copy_AAs 

int codon_2_num(char* codon)
{
	//fprintf(stderr, "codon_2_num needs fixing.\n");
	int num = 0;
	num += nuc_2_num(codon[0]) * 16 + 
			nuc_2_num(codon[1]) * 4 + 
			nuc_2_num(codon[2]);

	return(num);
}

void free_AAs_codons_mem(char** AAs_codons)
{
	for (int i_trip = 0; i_trip < 64; i_trip++)
	{
		delete[] AAs_codons[i_trip];
	} // i_trip loop.

	delete[] AAs_codons;
}

char* get_AA_code_per_codon(char* codon)
{
	//char** AAs = copy_AAs();

	int i_aa = codon_2_num(codon);

	 char* aa = new char[strlen(AAs[i_aa]) + 1];
	 strcpy(aa, AAs[i_aa]);

	 //free_AAs_codons_mem(AAs);

	return(aa);
}

char get_dna_pair(char nuc)
{
	if(toupper(nuc) == 'A')
	{
		return('T');
	}
	else if(toupper(nuc) == 'C')
	{
		return('G');
	}
	else if(toupper(nuc) == 'G')
	{
		return('C');
	}
	else if(toupper(nuc) == 'T')
	{
		return('A');
	}
	else
	{
		return('N');
	}
}

bool check_rna_pairing(char nuc1, char nuc2)
{
	if(toupper(nuc1) == 'A')
	{
		if(toupper(nuc2) == 'U' ||
			toupper(nuc2) == 'T')
		{
			return(true);
		}
	}

	if(toupper(nuc1) == 'C')
	{
		if(toupper(nuc2) == 'G')
		{
			return(true);
		}
	}

	if(toupper(nuc1) == 'G')
	{
		if(toupper(nuc2) == 'U' || 
			toupper(nuc2) == 'T' || 
			toupper(nuc2) == 'C')
		{
			return(true);
		}
	}

	if(toupper(nuc1) == 'T')
	{
		if(toupper(nuc2) == 'A' ||
			toupper(nuc2) == 'G')
		{
			return(true);
		}
	}

	if(toupper(nuc1) == 'U')
	{
		if(toupper(nuc2) == 'A' ||
			toupper(nuc2) == 'G')
		{
			return(true);
		}
	}

	return(false);
}

bool is_valid_nuc_char(char nuc)
{
	if(toupper(nuc) == 'A' ||
		toupper(nuc) == 'C' ||
		toupper(nuc) == 'G' ||
		toupper(nuc) == 'U' ||
		toupper(nuc) == 'T')
	{
		return(true);
	}
	else
	{
		return(false);
	}
}

int nuc_2_num(char nuc)
{
	if(toupper(nuc) == 'A')
	{
		return(0);
	}
	else if(toupper(nuc) == 'C')
	{
		return(1);
	}
	else if(toupper(nuc) == 'G')
	{
		return(2);
	}
	else if(toupper(nuc) == 'U' ||
		toupper(nuc) == 'T')
	{
		return(3);
	}
	else
	{
		return(4);
		//fprintf(stderr, "Cannot convert %c to number.\n", nuc);
		//exit(1);
	}
}

char num_2_nuc(int num)
{
	char nucs[] = "ACGT-";

	char nuc = nucs[num];

	return(nuc);
}

// Get the rna nucleotide corresponding to a given dna nucleotide.
char get_transcribed_rna_nuc_per_dna_nuc(char dna_nuc)
{
	if(toupper(dna_nuc) == 'A' || 
		toupper(dna_nuc) == 'C' || 
		toupper(dna_nuc) == 'G')
	{
		return(toupper(dna_nuc));
	}
	else if(toupper(dna_nuc) == 'T')
	{
		return('U');
	}
	else
	{
		return(toupper(dna_nuc));
	}
}

char get_transcribing_dna_nuc_per_rna_nuc(char rna_nuc)
{
	if(toupper(rna_nuc) == 'A' || 
		toupper(rna_nuc) == 'C' || 
		toupper(rna_nuc) == 'G')
	{
		return(toupper(rna_nuc));
	}
	else if(toupper(rna_nuc) == 'U')
	{
		return('T');
	}
	else
	{
		return(toupper(rna_nuc));
	}
}

int* get_per_char_as_nuc_2_num_coding_array()
{
	// code all non-nucleotide values as 4, others as a:0, c:1, g:2, t:3.
	int* per_nuc_val = new int[256];
	for(int i = 0; i < 256; i++)
	{
		per_nuc_val[i] = 4;
	} // i loop.
	per_nuc_val[(int)'A'] = 0;
	per_nuc_val[(int)'C'] = 1;
	per_nuc_val[(int)'G'] = 2;
	per_nuc_val[(int)'T'] = 3;
	per_nuc_val[(int)'U'] = 3;
	per_nuc_val[(int)'a'] = 0;
	per_nuc_val[(int)'c'] = 1;
	per_nuc_val[(int)'g'] = 2;
	per_nuc_val[(int)'t'] = 3;
	per_nuc_val[(int)'u'] = 3;

	return(per_nuc_val);
}

char get_complementary_dna_nuc_per_dna_nuc(char dna_nuc)
{
	if(toupper(dna_nuc) == 'A')
	{
		return('T');
	}
	else if(toupper(dna_nuc) == 'C')
	{
		return('G');
	}
	else if(toupper(dna_nuc) == 'G')
	{
		return('C');
	}
	else if(toupper(dna_nuc) == 'T' ||
		toupper(dna_nuc) == 'U')
	{
		return('A');
	}
	else
	{
		return('N');
	}
}

char get_complementary_rna_nuc_per_rna_nuc(char rna_nuc)
{
	if(toupper(rna_nuc) == 'A')
	{
		return('U');
	}
	else if(toupper(rna_nuc) == 'C')
	{
		return('G');
	}
	else if(toupper(rna_nuc) == 'G')
	{
		return('C');
	}
	else if(toupper(rna_nuc) == 'U')
	{
		return('A');
	}
	else
	{
		return('N');
	}
}

/*

A 
A adenine 

C 
C cytosine 

G 
G guanine 

T 
T thymine 

U 
U uracil 
   

R 
A or G purine 

Y 
C or T (U) pyrimidine 
   

M 
A or C amino 

K 
G or T (U) keto 

S 
C or G strong (3 H bonds) 

W 
A or T (U) weak (2 H bonds) 
   

B 
C or G or T (U) not A 

D 
A or G or T (U) not C 

H 
A or C or T (U) not G 

V 
A or C or G not T (U) 
   

N 
A or C or G or T (U) any nucleotide 

*/
void map_nuc_IUPAC_code(char raw_nuc, 
						char &trans_nuc, 
						int &num, 
						bool& force_unpaired)
{
	fprintf(stderr, "CODEBASE NUCLEOTIDE CODING IS NOT COMPATIBLE WITH IUPAC CODING.\n");
	exit(1);

	if(raw_nuc == 'a' || raw_nuc == 'c' || raw_nuc == 'g' || raw_nuc == 'u' || raw_nuc == 't')
	{
		force_unpaired = true;
	}
	else
	{
		force_unpaired = false;
	}

	if (toupper(raw_nuc) == 'A') 
	{
		trans_nuc = raw_nuc;
		num=1;
	}
	else if(toupper(raw_nuc) == 'B')
	{
		trans_nuc = 'N';
		num = 0;
	}
	else if(toupper(raw_nuc) == 'C')
	{
		trans_nuc = raw_nuc;
		num = 2;
	}
	else if(toupper(raw_nuc) == 'D')
	{
		trans_nuc = 'N';
		num = 0;
	}
	else if(toupper(raw_nuc) == 'G')
	{
		trans_nuc = raw_nuc;
		num = 3;
	}
	else if(toupper(raw_nuc) == 'H')
	{
		trans_nuc = 'N';
		num = 0;
	}
	else if(toupper(raw_nuc) == 'I')
	{
		trans_nuc = 'N';
		num = 0;
	}
	else if(toupper(raw_nuc) == 'K')
	{
		trans_nuc = 'N';
		num = 0;
	}
	else if(toupper(raw_nuc) == 'M')
	{
		trans_nuc = 'N';
		num = 0;
	}
	else if(toupper(raw_nuc) == 'N')
	{
		trans_nuc = 'N';
		num = 0;
	}
	else if(toupper(raw_nuc) == 'R')
	{
		trans_nuc = 'N';
		num = 0;
	}
	else if(toupper(raw_nuc) == 'S')
	{
		trans_nuc = 'N';
		num = 0;
	}
	else if(toupper(raw_nuc) == 'T')
	{
		trans_nuc = raw_nuc;
		num = 4;
	}
	else if(toupper(raw_nuc) == 'U')
	{
		trans_nuc = raw_nuc;
		num = 4;
	}
	else if(toupper(raw_nuc) == 'V')
	{
		trans_nuc = 'N';
		num = 0;
	}
	else if(toupper(raw_nuc) == 'W')
	{
		trans_nuc = 'N';
		num = 0;
	}
        else if(toupper(raw_nuc) == 'X')
        {
                trans_nuc = 'N';
                num = 0;
        }
	else if(toupper(raw_nuc) == 'Y')
	{
		trans_nuc = 'N';
		num = 0;
	}
	else
	{
		trans_nuc = 'N';
		num = 0;
	}

	if(num == 0)
	{
		//printf("Found %c\n", raw_nuc);
		//getc(stdin);
	}
}

void reverse_complement_seq(char* seq)
{
	int l_seq = strlen(seq);
	for(int i = 0; i < l_seq ; i++)
	{
		//char nuc_at_i = seq[i];
		seq[i] = get_complementary_dna_nuc_per_dna_nuc(seq[i]);
	} // i loop.

	char* reverted_seq = t_string::revert(seq);
	t_string::copy(seq, reverted_seq);
	delete [] reverted_seq;
}