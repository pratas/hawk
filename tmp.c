#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <ctype.h>

typedef unsigned long long ULL;
typedef unsigned char Char;
typedef unsigned char Base;
typedef unsigned char Symbol;

#define DEFAULT_HEADER_SIZE     512
#define DEFAULT_READ_SIZE       1024
#define READ_LEFT_GUARD         32

typedef struct{
  Char  *header1[2];
  Char  *bases;
  Char  *header2;
  Char  *scores;
  unsigned headerMaxSize;
  unsigned readMaxSize;
  unsigned solidData;
  unsigned header2Present;
  unsigned skipNs;
  uint8_t  lowestScore;
  }
Read;

Read *CreateRead(unsigned headerMaxSize, unsigned readMaxSize)

        {
        Read *read;

        if(!(read = (Read *)Calloc(1, sizeof(Read))))
                {
                fprintf(stderr, "Error: out of memory\n");
                exit(1);
                }

        /* The header1 from the previous read */
        if(!(read->header1[0] = (Char *)Calloc(headerMaxSize + READ_LEFT_GUARD,
          sizeof(Char))))
                {
                fprintf(stderr, "Error: out of memory\n");
                exit(1);
                }

        /* The header1 from the current read */
        if(!(read->header1[1] = (Char *)Calloc(headerMaxSize + READ_LEFT_GUARD,
          sizeof(Char))))
                {
                fprintf(stderr, "Error: out of memory\n");
                exit(1);
                }

        if(!(read->header2 = (Char *)Calloc(headerMaxSize + READ_LEFT_GUARD,
          sizeof(Char))))
                {
                fprintf(stderr, "Error: out of memory\n");
                exit(1);
                }

        read->headerMaxSize = headerMaxSize;

        if(!(read->bases = (Char *)Calloc(readMaxSize + READ_LEFT_GUARD,
          sizeof(Char))))
                {
                fprintf(stderr, "Error: out of memory\n");
                exit(1);
                }

        if(!(read->scores = (Char *)Calloc(readMaxSize + READ_LEFT_GUARD,
          sizeof(Char))))
                {
                fprintf(stderr, "Error: out of memory\n");
                exit(1);
                }

        read->readMaxSize = readMaxSize;

        read->header1[0] += READ_LEFT_GUARD;
        read->header1[1] += READ_LEFT_GUARD;
        read->header2 += READ_LEFT_GUARD;
        read->bases += READ_LEFT_GUARD;
        read->scores += READ_LEFT_GUARD;

        read->solidData = 0;
        read->header2Present = 0;
        read->skipNs = 0;
        read->lowestScore = (Char)255;

        return read;
        }


Read *GetRead(FILE *fp, Read *read)
        {
        int n, c = fgetc(fp);

        if(c == EOF)
                return NULL;

        /* Check if the initial '@' is present in the header line */
        if(c != '@')
                {
                fprintf(stderr, "Error: failed to get the initial '@' character\n");
                exit(1);
                }

        if(!fgets((char *)read->header1[1], read->headerMaxSize, fp))
                {
                fprintf(stderr, "Error: unexpected end of file\n");
                exit(1);
                }

        if(!fgets((char *)read->bases, read->readMaxSize, fp))
                {
                fprintf(stderr, "Error: unexpected end of file\n");
                exit(1);
                }

        if(!fgets((char *)read->header2, read->headerMaxSize, fp))
                {
                fprintf(stderr, "Error: unexpected end of file\n");
                exit(1);
                }

        if(!fgets((char *)read->scores, read->readMaxSize, fp))
                {
                fprintf(stderr, "Error: unexpected end of file\n");
                exit(1);
                }

        if(read->solidData)
                {
                n = 1;
                while(read->bases[n] != '\n')
                        {
                        read->bases[n] = PseudoCsToBase(read->bases[n - 1], read->bases[n]);
                        n++;
                        }

                }

        return read;
        }



int main(int argc, char *argv[]){
  FILE *F=fopen("B.fastq", "r");
  Read R = CreateRead(DEFAULT_HEADER_SIZE, DEFAULT_READ_SIZE); 
  int long long nReads = 0, fileSize = 0;
      
  while(GetRead(F, R)){
    fileSize += strlen((char *)read->header1[1]) + 1;
    fileSize += strlen((char *)read->header2);
    fileSize += strlen((char *)read->bases);
    fileSize += strlen((char *)read->scores);
    }

  fprintf("filesize: %llu\n", fileSize);
  return 0;
  }
