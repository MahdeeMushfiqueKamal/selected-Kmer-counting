
#include <iostream>
#include <fstream>

using namespace std;



uint64_t cal(string str_k_mer)
{
    const char* char_array = str_k_mer.c_str();
    uint64_t k_mer = 0;
    for(int j=0; char_array[j]!='\0'; j++)
    {
        switch(char_array[j])
        {
            case 'A':
                    k_mer = k_mer<<2; // A=00
                    break;
            case 'T':
                    k_mer = k_mer<<2;
                    k_mer = k_mer | 1;   //T=01
                    break;
            case 'G':
                    k_mer = k_mer<<2;  //G=10
                    k_mer = k_mer | 2;
                    break;
            case 'C':
                    k_mer = k_mer<<2; //C=11
                    k_mer = k_mer | 3;
                    break;
        }
    }
    return k_mer;
}




int64_t addChar(int64_t k_mer, char ch)
{
    switch(ch)
    {
        case 'A':
                k_mer = k_mer<<2; // A=00
                break;
        case 'T':
                k_mer = k_mer<<2;
                k_mer = k_mer | 1;   //T=01
                break;
        case 'G':
                k_mer = k_mer<<2;  //G=10
                k_mer = k_mer | 2;
                break;
        case 'C':
                k_mer = k_mer<<2; //C=11
                k_mer = k_mer | 3;
                break;
    }
    return k_mer;

}


uint64_t reverse_compliment(uint64_t cal_kmer, int kmer_length)
{
    uint64_t k_mer = 0;
    uint64_t mask = 3;

    for(int i=0; i<kmer_length; i++)
    {
        switch(cal_kmer & mask)
        {

        case 0:
            k_mer = k_mer<<2;
            k_mer = k_mer | 1;   //A=00->T=01
            break;
        case 1:
            k_mer = k_mer<<2; //T=01->A=00
            break;
        case 2:
            k_mer = k_mer<<2; //G=10->C=11
            k_mer = k_mer | 3;
            break;
        case 3:
            k_mer = k_mer<<2;  //C=11->G=10
            k_mer = k_mer | 2;
            break;

        }
        cal_kmer=cal_kmer>>2;
    }
    return k_mer;
}


void read_kmer(bool canonicalize) {

    std::ifstream fasta_file("lhg22L20MC5x.fa");
    int64_t currentKmer = 0;
    int current_len = 0;
    int chunk_size = 500;
    char arr_chunk[chunk_size];
    int kmerLen = 22;
    int occurance = 0;

    bool isInHeader = false;





    char sample1[]="TTCAAAGCTGCTCTGTAAAAAG";
    char sample2[]="CTTTTTACAGAGCAGCTTTGAA";
//    char sample3[]="ATGAATGGAATCGAATGGAATC";
    int64_t sample_kmer1 = cal(sample1);
    int64_t sample_kmer2 = cal(sample2);
//    int64_t sample_kmer3 = cal(sample3);


    int64_t mask_1, mask_2;
        if(2*kmerLen>30){
            mask_1 = (1 << 30);
            mask_1 = ~(mask_1 << (2*(kmerLen-1)-30));
            mask_2 = (1 << 30);
            mask_2 = ~(mask_2 << ((2*(kmerLen-1)-30)+1));
            //cout<<"inside if"<<endl;

        }
        else{
            mask_1 = ~(1<<(2*(kmerLen-1)));
            mask_2 = ~(1<<((2*(kmerLen-1))+1));
            //cout<<"inside else"<<endl;
        }

    int chunk_count = 0;
    int kmer_count = 0;












    while(!fasta_file.eof())
    {
        chunk_count++;
        fasta_file.read(arr_chunk, chunk_size);
        int i = 0;

        while(i<chunk_size)
        {
            if(arr_chunk[i]==EOF)
                break;

            char ch = arr_chunk[i];
            if(ch == '>')
            {
                isInHeader = true;
                currentKmer = 0;
                current_len = 0;
            }
            else if (isInHeader==true && ch=='\n')
            {
                isInHeader = false;
                currentKmer = 0;
                current_len = 0;
            }

            if(isInHeader)
            {
                i+=1;
                continue;
            }

            if(ch=='\n' || ch=='\r' || ch==' ')
            {
                i+=1;
                continue;
            }
            if(ch=='N')
            {
                currentKmer = 0;
                current_len = 0;
                i+=1;
                continue;
            }
            else
            {
                kmer_count++;
                if(current_len < kmerLen)
                {
                    currentKmer = addChar(currentKmer, ch);
                    current_len++;
                }
                else
                {
                    currentKmer = currentKmer & mask_1;
                    currentKmer = currentKmer & mask_2;
                    currentKmer = addChar(currentKmer, ch);
                    //current_len++;
                }
                //cout<<currentKmer<<endl;

                if(current_len==22)
                {
                    //GOT KMER
                    if(canonicalize) reverse_compliment(currentKmer, kmerLen);
                }

                if(currentKmer==sample_kmer1 && current_len==22 || currentKmer==sample_kmer2 && current_len==22)
                {
                    //kmer_count++;
                    occurance+=1;

                }


            }




            i+=1;


        }

    }
    cout<<occurance<<endl;
    //cout<<chunk_count<<endl;


}


int main() {

    //cout<<cal("AATGGAATCAAATGGAATCATC")<<endl;
    //cout<<cal("ATGAATGGAATCGAATGGAATC")<<endl;
    read_kmer(true);


    return 0;
}
