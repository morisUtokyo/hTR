// Scan minimap2 cigar file

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include <getopt.h>

#define PRINT_nearest_SNV

#define BLK 4096
#define DISTANCE 2000
#define MIN_CS_TH 990
#define MINIMUM_FREQ 4
// We retain such significant single nucleotide variants with k occurrences in n reads that the probability of observing k or more occurrences is smaller than 5% level of significance when a substitution is assumed to be observed at random with probability p (e.g., 0.05%). For example, when n=3000 and p=0.05%, we can set k=4.

struct hap{
    char type[BLK];
    int  freq;
    char ID[10];
};
#define hapTblSize 100
struct hap hapTbl[hapTblSize];
int hapTblIndex = 0;

int add2hapTbl(char *oneHap, char *tmpID){
    for(int i=0; i<hapTblIndex; i++){
        if(strcmp(hapTbl[i].type, oneHap) == 0){    // oneHap has been already added
            if(strcmp(hapTbl[i].ID, tmpID) != 0){   // oneHap occurs first for tmpID
                hapTbl[i].freq++;
                strcpy(hapTbl[i].ID, tmpID);
                return 1;
            }else{ // oneHap for tmpID has been found
                return 0;
            }
        }
    }
    // absence of one Hap
    strcpy( hapTbl[hapTblIndex].type, oneHap);
    hapTbl[hapTblIndex].freq = 1;
    strcpy( hapTbl[hapTblIndex].ID,   tmpID);
    hapTblIndex++;
    return 1;
}
void dump_hapTbl(){
    //printf("\nhapTblIndex=%d\n", hapTblIndex);
    for(int i=0; i<hapTblIndex; i++)
        printf("%d\t%s\n", hapTbl[i].freq, hapTbl[i].type);
}


int** SNV_Tbl;   // (loc, freq)
int SNV_TblIndex;
#define SNV_TblSize 10000

void init_SNV_Tbl(){
    SNV_Tbl = (int**)malloc(sizeof(int*) * SNV_TblSize);   // (loc, freq)
    for(int i=0; i<SNV_TblSize; i++){
        SNV_Tbl[i] = (int*)malloc(sizeof(int) * 2);
        if(SNV_Tbl[i] == NULL) exit(EXIT_FAILURE);
    }
    SNV_TblIndex = 0;
}

void free_SNV_Tbl(){
    for(int i=0; i<SNV_TblSize; i++) free(SNV_Tbl[i]);
    free(SNV_Tbl);
}

void sort_SNV_Tbl() {
    for(int i = 1; i < SNV_TblIndex; i++) {
        if(SNV_Tbl[i-1][1] < SNV_Tbl[i][1]) { // Compare the frequencies
            int j = i;
            int tmp_pos  = SNV_Tbl[i][0];
            int tmp_freq = SNV_Tbl[i][1];
            do {
                SNV_Tbl[j][0] = SNV_Tbl[j-1][0];
                SNV_Tbl[j][1] = SNV_Tbl[j-1][1];
                j--;
            } while (j > 0 && SNV_Tbl[j-1][1] < tmp_freq);
            SNV_Tbl[j][0] = tmp_pos;
            SNV_Tbl[j][1] = tmp_freq;
        }
    }
}

void add2SNV_Tbl(int readPos){
    for(int i=0; i<SNV_TblIndex; i++){
        if(SNV_Tbl[i][0] == readPos){
            SNV_Tbl[i][1]++;
            return;
        }
    }
    SNV_Tbl[SNV_TblIndex][0] = readPos;
    SNV_Tbl[SNV_TblIndex][1] = 1;
    SNV_TblIndex++;
    if(SNV_TblIndex%10 == 0) sort_SNV_Tbl();
}

int freq(int readPos){
    for(int i=0; i<SNV_TblIndex; i++){
        if(SNV_Tbl[i][0] == readPos){
            return(SNV_Tbl[i][1]);
        }
    }
    return(0);
}

void dump_SNV_Tbl(){
    printf("SNV_TblIndex=%d\n", SNV_TblIndex);
    for(int i=0; i<SNV_TblIndex; i++)
        printf("%d\t%d\n", SNV_Tbl[i][0], SNV_Tbl[i][1]);
}

int length_cs(char *cs){
    // Compute the length of a cs SAM/PAF tag
    // Example    :6-ata:10+gtc:4*at:3
    // CGATCCATAAATAGAGTAG---GAAAGCA
    // ||||||   ||||||||||   ||| |||
    // CGATCC---AATAGAGTAGGTCGAATGCA
    // : #matches, - deletion, + insertion, * one mismatch
    int readPos = 0;
    char numString[BLK];
    int j = 0; // position in numString
    
    for(int i=0; ; i++){
        char prev_c;
        char c=cs[i];
        if(c==':' || c=='-' || c=='+' || c=='*' || c=='\0'){
            if(j > 0){
                numString[j]='\0';
                int k;
                switch(prev_c){
                    case ':':
                        sscanf(numString, "%d", &k); readPos += k; break;
                    case '-':   readPos += j; break;
                    case '+':   break;
                    case '*':   readPos += 1; break;
                }
            }
            if(c=='\0') break;  // The end of the cs string
            else{   j = 0;  prev_c = c; }
        }else
            numString[j++] = c;
    }
    return readPos;
}


int parse_cs(char *pre_post, int refPos, char *cs, int initial, char *haplotype){
    sprintf(haplotype, "");
    // Parse a cs SMD/PAF tag
    
    int pre;
    if(strcmp(pre_post, "pre")==0)  pre = 1;
    if(strcmp(pre_post, "post")==0) pre = 0;
    char pre_SNV[BLK], post_SNV[BLK];
    strcpy(pre_SNV, ""); strcpy(post_SNV, "");
    int first_SNV_flag = 1;

    if(length_cs(cs) < MIN_CS_TH)   return 0;
    
    int readPos = refPos;
    char numString[BLK];
    int j = 0; // position in numString
    
    for(int i=0; ; i++){
        char prev_c;
        char c=cs[i];
        if(c==':' || c=='-' || c=='+' || c=='*' || c=='\0'){
            if(j > 0){
                numString[j]='\0';
                int k;
                switch(prev_c){
                    case ':':   // matches
                        sscanf(numString, "%d", &k);
                        readPos += k;
                        break;
                    case '-':   // deletion
                        readPos += j;
                        break;
                    case '+':   // insertion
                        break;
                    case '*':   // one mismatch SNV
                        readPos += 1;
                        if(initial == 1){
                            add2SNV_Tbl(readPos);
                        }else{
                            int frequency = freq(readPos);
                            if(MINIMUM_FREQ <= frequency){
                                sprintf(haplotype, "%s%d", haplotype, readPos);
                                // Set pre_SNV to the last SNV
                                if(pre == 1)
                                    sprintf(pre_SNV, "%d", readPos);
                                if(pre == 0 && first_SNV_flag == 1){
                                    // Assign the first SNV to post_SNV
                                    sprintf(post_SNV, "%d", readPos);
                                    first_SNV_flag = 0;
                                }
                            }
                        }
                        break;
                }
            }
            if(c=='\0') break;  // The end of the cs string
            else{
                j = 0;
                prev_c = c;
            }
        }else{
            numString[j++] = c;
        }
    }
#ifdef PRINT_nearest_SNV
    if(pre == 1)
        sprintf(haplotype, "%s", pre_SNV);
    else
        sprintf(haplotype, "%s", post_SNV);
#endif
    return 1;
}

int check_position(char *pre_post, int chrNum, int chrID, int beginTR, int startPos, int endTR){
    
    int flag=0;
    if(strcmp(pre_post, "pre")==0){
        // The previous sequence alignment must begin after beginTR - DISTANCE
        if(chrNum == chrID &&
           beginTR - DISTANCE <= startPos && startPos <= beginTR)
            flag = 1;
    }else if(strcmp(pre_post, "post")==0){
        // The post sequence alignment must begin after endTR
        if(chrNum == chrID &&
           endTR <= startPos && startPos <= endTR + DISTANCE)
            flag = 1;
    }
    return flag;
}

int scan_cigar_cs(char *others, char *cigar, char *cs){
    char head[BLK], tail[BLK];
    
    int fed_cigar = 0;
    int fed_cs = 0;
    int last = 0;

    while(strcmp(others, "") != 0){
        if(sscanf(others, "%s\t%[^\n]", head, tail) < 2 ){
            // return 2 if others match the pattern
            strcpy(head, others);
            last = 1;
        }
        //printf("head=%s tail=%s\n", head, tail);
        
        if(fed_cigar == 0){
            if(sscanf(head, "cg:Z:%s", cigar) == 1)
                fed_cigar = 1;
        }
        if(fed_cs == 0){
            if(sscanf(head, "cs:Z:%s", cs) == 1)
                fed_cs = 1;
        }
        if(last == 0)
            strcpy(others, tail);
        else
            break;
    }
    if(fed_cigar == 1 && fed_cs == 1){
        //printf("cigar=%s cs=%s\n", cigar, cs);
        return(1);
    }else
        return(0);
}

int main(int argc, char *argv[])
{    
    char inputFile[BLK];    // For the input file name
    int inputFile_given = 0;
    int opt;
    int chrNum, beginTR, endTR;
    while ((opt = getopt(argc, argv, "f:c:b:e:")) != -1) {
        switch(opt){
            case 'f':
                strcpy(inputFile,optarg);
                inputFile_given = 1;
                break;
            case 'c':
                sscanf(optarg, "%d", &chrNum);
                //fprintf(stderr, "%d\n", chrNum);
                break;
            case 'b':
                sscanf(optarg, "%d", &beginTR);
                //fprintf(stderr, "%d\n", beginTR);
                break;
            case 'e':
                sscanf(optarg, "%d", &endTR);
                //fprintf(stderr, "%d\n", endTR);
                break;
            default:
                fprintf(stderr, "Usage: hap -f <minimap2 CIGAR file> -c <chr number> -b <begin of TR> -e <end of TR> \n");
                exit(EXIT_FAILURE);
        }
    }
    if(inputFile_given == 0){
        fprintf(stderr, "Input file is not given.\n");
        exit(EXIT_FAILURE);
    }
 
    init_SNV_Tbl();
    char *s = (char *)malloc(sizeof(char)*BLK);
    char ID[BLK], pre_post[10], strand[10], cigar[BLK], cs[BLK], haplotype[BLK], others[BLK];
    int  readID, chrID, startPos, d1, d2, d3, d4;
    
    // Count the frequency of each SNV
    FILE *fp = fopen(inputFile, "r");
    while (fgets(s, BLK, fp) != NULL) {
        //sscanf(s, "%[^,],read%d,%s %s chr%d %d cg:Z:%s cs:Z:%s", ID, &readID, pre_post, strand, &chrID, &startPos, cigar, cs);
        // 1st ID, readID, pre_post
        // 2nd, 3rd, 4th are skipped by using dummy variables
        // 5th strand +/-
        // 6th chr1 etc.
        // 7th skipped by %*c
        // 8th startPos
        sscanf(s, "%[^,],read%d,%s\t%d\t%d\t%d\t%s\tchr%d\t%d\t%d\t%[^\n]", ID, &readID, pre_post, &d1, &d2, &d3, strand, &chrID, &d4, &startPos, others);
        // parse cigar and cs from other features
        int fed_cigar_cs = scan_cigar_cs(others, cigar, cs);
        //printf("ID=%s readID=%d %s chr%d startPos=%d cigar=%s cs=%s\n", ID, readID, pre_post, chrID, startPos, cigar, cs);
        
        if(fed_cigar_cs == 1){
            int flag = check_position(pre_post, chrNum, chrID, beginTR, startPos, endTR);
            if(flag == 1){
                parse_cs(pre_post, startPos, cs, 1, haplotype);  // 1 = Initialize SNV_Tbl
            }
        }
    }
    sort_SNV_Tbl();
    //dump_SNV_Tbl();
    fclose(fp);
    
    // Generate haplotypes
    fp = fopen(inputFile, "r");
    int numQualified = 0;
    int readIDprev = -1;
    char haplotypeAll[BLK];
    char prevID[BLK];
    while (fgets(s, BLK, fp) != NULL) {
        //sscanf(s, "%[^,],read%d,%s %s chr%d %d cg:Z:%s cs:Z:%s", ID, &readID, pre_post, strand, &chrID, &startPos, cigar, cs);
        // 1st ID, readID, pre_post
        // 2nd, 3rd, 4th are skipped by using dummy variables
        // 5th strand +/-
        // 6th chr1 etc.
        // 7th skipped by %*c
        // 8th startPos
        sscanf(s, "%[^,],read%d,%s\t%d\t%d\t%d\t%s\tchr%d\t%d\t%d\t%[^\n]", ID, &readID, pre_post, &d1, &d2, &d3, strand, &chrID, &d4, &startPos, others);
        // parse cigar and cs from other features
        int fed_cigar_cs = scan_cigar_cs(others, cigar, cs);

        if(readIDprev != readID){
            if(numQualified == 2){  // Print all reads
                printf("%s\t%d\t%s\n", prevID, readIDprev, haplotypeAll);
            }
            sprintf(haplotypeAll, "");
            readIDprev = readID;
            numQualified = 0;
        }
        strcpy(prevID, ID);
        int flag = check_position(pre_post, chrNum, chrID, beginTR, startPos, endTR);   // The position is consistent with the given range.
        if(flag == 1){
            int qualified = parse_cs(pre_post, startPos, cs, 0, haplotype);
            if(qualified == 1){ // The surrounding string is long enough (~1Kb).
                numQualified++;
                if(strcmp(pre_post, "pre")==0)
                    sprintf(haplotypeAll, "%s|", haplotype);
                if(strcmp(pre_post, "post")==0)
                    sprintf(haplotypeAll, "%s%s", haplotypeAll, haplotype);
            }
        }
        
    }
    dump_hapTbl();
    fclose(fp);
    
    //dump_SNV_Tbl();
    free_SNV_Tbl();
    free(s);

    return EXIT_SUCCESS;
}
