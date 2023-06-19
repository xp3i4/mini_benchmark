#include <seqan/seq_io.h>
#include "shape_extend.h"
using namespace seqan; 


int getHashStatistics(std::vector<int> & v1, double density, double time, std::vector<double> & res, LShape shape, CharString cs, CharString cs2)
{

    std::vector<Pair<int, int> > bin_size_count;
    std::vector<int> tmp;
    uint64_t pre = v1[0];
    uint64_t bin_size = 0;
    uint64_t all = 0;
    double q1 = 0, q2 = 0, q3 = 0, qmax = 0, qmax_v=0, q4 = 0;
    std::vector<int> v;
    for (int i = 0; i < v1.size(); i++)
    {
        if (pre != v1[i])     
        {
            v.push_back(pre);
        }
        pre = v1[i];
    }
    std::sort(v.begin(), v.end());
    for (unsigned i = 0; i < v.size(); i++)
    {
        //dout << v[i] << pre << "\n";
        if (v[i] != pre)
        {
            tmp.push_back(bin_size);
            pre = v[i];
            bin_size = 0;
        }
        if (qmax < bin_size)
        {
            qmax = bin_size;
            qmax_v = pre;
        }
        pre = v[i];
        bin_size++;
    }
    std::sort(tmp.begin(), tmp.end());
    pre = tmp[0];
    int count = 0;
    for (unsigned i = 0; i < tmp.size(); i++)
    {
        if(tmp[i] != pre)
        {
            bin_size_count.push_back(Pair<int,int>(pre, count));
            all += count;
            count = 0;
        }
        pre = tmp[i];
        count++;
    }
    bin_size_count.push_back(Pair<int,int>(pre, count));
    all += count;

    uint64_t s = 0;
    double count_0 = 0;
    for (double i = 1; i < bin_size_count.size(); i++) 
    {
        if (s < all * 0.25 && s + bin_size_count[i].i2 >= all * 0.25)
        {
            q1 = bin_size_count[i].i1;
        }
        if (s < all * 0.5 && s + bin_size_count[i].i2 >= all * 0.5)
        {
            q2 = bin_size_count[i].i1;
        }
        if (s < all * 0.75 && s + bin_size_count[i].i2 >= all * 0.75) 
        {
            q3 = bin_size_count[i].i1;
        }
        if (s < all * 0.95 && s + bin_size_count[i].i2 >= all * 0.95) 
        {
            q4 = bin_size_count[i].i1;
        }
        s += bin_size_count[i].i2;
        //std::cout << bin_size_count[i].i1 * double(1000000) / v.size() <<","<< bin_size_count[i].i2 /float(all)<<","<< shape.span <<","<< shape.weight << "," << cs << ","<< cs2 <<"\n";
    }
    //dout << "all= " << all << "\n";
    double p_u = double(1) / (1 << shape.weight * 2); //Probability of uniform distribution
    double d_kl = 0; //KL divergence
    for (double i = 0; i < tmp.size(); i++) 
    {
        double p_x = double(tmp[i]) / v.size();
        if (p_x != 0)
        {
            //dout << "p_x" << p_x << p_u << "\n";
            d_kl += p_x * std::log(p_x / p_u);
        }

    }
    res.push_back(q1 * 1000000 / v.size());
    res.push_back(q2 * 1000000 / v.size());
    res.push_back(q3 * 1000000 / v.size());
    res.push_back(q4 * 1000000 / v.size());
    res.push_back(density);
    res.push_back(d_kl);
    res.push_back(time);

    //std::cout << "all= "  << v.size() << " " << qmax << " " << all << " " << qmax_v << "\n";
    //dout << "Quantiles " << q1 * 1000000 / v.size() << q2 * 1000000 / v.size() << q3 * 1000000 / v.size() << qmax * 1000000 / v.size() <<  d_kl << "\n";
    return 0;
}

int main(int argc, char const ** argv)
{
    double time = sysTime();
    std::cerr << std::fixed << std::setprecision(2);
    SeqFileIn rFile;
    StringSet<CharString> read_id;
    StringSet<String<seqan::Dna5> > reads;
    open(rFile, toCString("t.fa"));
    readRecords(read_id, reads, rFile, 500);
    LShape shape(0);
    shape.span = 25;
    shape.weight = 10;
    std::vector<int> v1, v2;
    std::vector<std::pair<unsigned, unsigned> > shape_weight_pairs;
    shape_weight_pairs.push_back(std::pair<unsigned, unsigned>(15,4));
    shape_weight_pairs.push_back(std::pair<unsigned, unsigned>(15,8));
    shape_weight_pairs.push_back(std::pair<unsigned, unsigned>(15,12));
    shape_weight_pairs.push_back(std::pair<unsigned, unsigned>(25,4));
    shape_weight_pairs.push_back(std::pair<unsigned, unsigned>(25,8));
    shape_weight_pairs.push_back(std::pair<unsigned, unsigned>(25,12));
    std::cout << "x,y,w,k,Type,Scheme" << "\n";
    for (unsigned s = 0; s < length(shape_weight_pairs); s++)
    {
        double t = 0, t1 = 0, t2 = 0;
        shape.span = shape_weight_pairs[s].first;
        shape.weight = shape_weight_pairs[s].second;
        double d1 = 0, d2 = 0;
        uint64_t pre1 = 0, pre2 = 0, ck = 0;
        for (unsigned i = 0; i < length(reads); i++)
        {
            hashInit(shape, begin(reads[i]));
            for (int j = 0; j < int(length(reads[i]) - shape.span - 1); j++)
            {
                hashNexth(shape, begin(reads[i]) + j);
                t = sysTime();
                hashNextX(shape, begin(reads[i]) + j);
                t1 += sysTime() - t;
                //std::cout << "x " << shape.XValue << " " << shape.x << "\n";
                v1.push_back(shape.XValue);
                if (pre1 != shape.XValue)
                {
                    d1++;
                    pre1 = shape.XValue;
                }
                t = sysTime();
                hashNextS(shape, begin(reads[i]) + j);
                t2 += sysTime() - t;
                v2.push_back(shape.XValue);
                if (pre2 != shape.XValue)
                {
                    d2++;
                    pre2 = shape.XValue;
                }
                ck++;
            }
        }
        d1 /= ck;
        d2 /= ck;
        std::vector<double> res1, res2;
        getHashStatistics(v1, d1, t1, res1, shape, "Rfd", "Lexico");
        //getHashStatistics(v1, d1, t1, res1, shape, "Rfd", "Random");
        getHashStatistics(v2, d2, t2, res2, shape, "Std", "Lexico");
        //getHashStatistics(v2, d2, t2, res2, shape, "Std", "Random");
        
        std::cout << std::fixed;
        std::cout << std::setprecision(2);
        std::cout << " & " << shape.span << ", "<< shape.weight << " && ";
        for (unsigned i = 0; i < length(res1); i++)
        {
            if (res1[i] < res2[i])
            {
                std::cout << res2[i] << " & \\bf{" << res1[i] << "}";
            }
            else 
            {
                std::cout << "\\bf{" << res2[i] << "} & " << res1[i];
            }
            if (i != length(res1) - 1)
            {
                std::cout << " && ";
            }
        }
        std::cout << "\\" << "\\ \n";
        
        clear(res1);
        clear(res2);
        clear(v1);
        clear(v2);
    }
    close(rFile);
    return 0;
}
