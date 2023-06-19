#include "shape_extend.h"

using namespace seqan;

//WARN!!:: Only odd Shape len is allowed if call the hash due to the double strand hash value
//All even shape len will be converted to len + 1
typedef Dna5 ShapeType;
std::vector<uint64_t> randoms;
void LShape::init_shape_parm (unsigned shape_span)
{
    span = shape_span ;
    weight = span - 8;
}
LShape::LShape(unsigned shape_span):
        hValue(0),
        crhValue(0),
        XValue(0),
        YValue(0),
        strand(0),
        leftChar(0),
        x(0)
{
    init_shape_parm(shape_span);
}

uint64_t randomInt64(uint64_t key, uint64_t mask)
{
    key = (~key + (key << 21)) & mask; 
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask; 
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask; 
    key = key ^ key >> 28;
    key = (key + (key << 31)) & mask;

    return key;
}

int randomOrder(std::vector<uint64_t> & s, unsigned len)
{
    if (s.size() != len)
    {
        s.resize(len);
        for (uint64_t i = 0; i < len; i++)
        {
            s[i] = i;
        }
        std::random_shuffle(s.begin(), s.end());

        for (uint64_t i = 0; i < len; i++)
        {
            //std::cout << "random " << s[i] << "\n";
        }
    }

    return 0;
}

void resize(LShape & me, unsigned new_span, unsigned new_weight)
{   
    me.span = new_span;
    me.weight = new_weight;
}   

static const uint64_t COMP4 = 3;
static const int  ordC = 3;
static const int pd[4]={1,1,-1,-1};

uint64_t getMask(unsigned bit)
{
    return (1ULL << bit) - 1;
}

uint64_t hashInit(LShape & me, TIterS it)
{
    SEQAN_ASSERT_GT((unsigned)me.span, 0u);
    randomOrder(randoms, 1 << (2 * me.weight));
    me.leftChar = 0;
    me.hValue = 0;
    me.crhValue = 0;
    me.leftChar = 0;
    me.x = me.leftChar - ordC;
    uint64_t k =0, count = 0; //COMP for complemnet value;
    while (count < me.span)
    {
        if (ordValue((ShapeType)*(it + k + count)) == 4)
        {
            k += count + 1;
            count = 0;
        }
        else
            count++;
    }
    unsigned bit = 2;
    for (unsigned i = 0; i < me.span - 1; ++i)
    {
        uint64_t val = ordValue ((ShapeType)*(it + k + i));
//        me.x += (int(val) << 1) - ordC;
        me.x += pd[val];
        me.hValue = (me.hValue << 2) + val;
        me.crhValue += ((COMP4 - val) << bit);
        bit += 2;
    }
    return k;
}
/*
 * this hashNext function is for index only collect mini hash value [minindex]
 * calculate hValue;
 */ 
uint64_t hashNexth(LShape & me, TIterS it)
{
    SEQAN_ASSERT_GT((unsigned)me.span, 0u);
    uint64_t mask = getMask((me.span << 1) - 2);
    int v2 = ordValue((ShapeType)*(it + me.span - 1 ));
    me.hValue = ((me.hValue & mask) << 2) + v2;
    me.crhValue = ((me.crhValue >> 2) & mask) + 
                  ((COMP4 - v2) << ((me.span << 1) - 2));
    me.x += pd[v2] - pd[me.leftChar];
    //std::cout << v2 << " " << pd[v2] << " " << me.leftChar << " " << pd[me.leftChar] << " " << me.x  << "\n";
//    me.x += (v2 - me.leftChar) << 1;
    me.leftChar = ordValue((ShapeType)*(it));
    return me.x < 0 ? me.hValue : me.crhValue;
}
uint64_t hashNextX(LShape & me, TIterS it)
{
    SEQAN_ASSERT_GT((unsigned)me.span, 0u);
    uint64_t v1;
    unsigned span = me.span << 1, weight = me.weight << 1;
    uint64_t v2 = me.x < 0 ? me.hValue : me.crhValue;
    me.XValue = getMask(me.span << 1);
    uint64_t mask = (1ULL<<2 * me.weight) - 1;
    for (unsigned k = 64 - span; k <= 64 - weight; k += 2)
    {
        v1 = v2 << k >> (64 - weight);
        v1 = randomInt64(v1, mask);
        if(me.XValue > v1)
        {
            me.XValue = v1;
        }
    } 
    (void)it;
    return me.XValue; 
}
uint64_t hashNextS(LShape & me, TIterS it)
{
    SEQAN_ASSERT_GT((unsigned)me.span, 0u);
    uint64_t v1;
    unsigned span = me.span << 1, weight = me.weight << 1;
    uint64_t t = 0, v2;
    me.XValue = getMask(me.span << 1);
    v2 = me.hValue;
    uint64_t mask = (1ULL<<2 * me.weight) - 1;

    for (unsigned k = 64 - span; k <= 64 - weight; k += 2)
    {
        v1 = v2 << k >> (64 - weight);
        v1 = randomInt64(v1, mask);
        if(me.XValue > v1)
        {
            me.XValue = v1;
        }
    } 
    v2 = me.crhValue;
    for (unsigned k = 64 - span; k <= 64 - weight; k += 2)
    {
        v1 = v2 << k >> (64 - weight);
        v1 = randomInt64(v1, mask);
        if(me.XValue > v1)
        {
            me.XValue = v1;
        }
    } 

    (void)it;
    return me.XValue; 
}


