#ifndef ONTOLOGYRULES_H
#define ONTOLOGYRULES_H

#include "../../ontology.h"
#include "../../node.h"

struct bottomFeature
{
    int IDind;  //my index in an array
    int posFeature;
    Node *noderef;
    Node *termRefTopDown;
    Node *termRefBottomUp;
    boost::dynamic_bitset<>  allChilds; //represents ids to more specific nodes
    boost::dynamic_bitset<>  allParents; //represents ids to more general nodes
    boost::dynamic_bitset<> exampleCovered;
    boost::dynamic_bitset<>  allSpecific;   //all specifics nodes
    boost::dynamic_bitset<>  allGeneral;    //all general nodes
    int level;  //level of scpecificity, i.e. root is 0, ...
    double actScore;
    double pathBestScore;
};

//represents new rule
struct newComplex
{
    double score;
    double pathScore;
    int coverPos;
    int coverNeg;
    double sig;
    std::vector<bottomFeature*> rules;
    //std::vector<bottomFeature*> candidates;
    bool operator<(const newComplex& rhs) const
    {
        return score > rhs.score;
    }
};

//represents new rule
struct newComplexStat
{
    double score;
    double pathScore;
    int coverPos;
    int coverNeg;
    double sig;
    std::vector<bottomFeature> rules;
};

class OntologyRules
{
public:
    OntologyRules(Ontology *ontoref, int ontoIndex, std::vector<Node*> *roots, std::vector<std::vector<boost::dynamic_bitset<> *> > *POSexamples, std::vector<std::vector<boost::dynamic_bitset<> *> > *NEGexamples);
    ~OntologyRules();
private:
    int ontoIndex;  //given by refOntologies
    Ontology *ontoref;
    std::vector<int> rootsIds;
    std::vector<bottomFeature> bottomRules;

    void initBottomRules(std::vector<std::vector<boost::dynamic_bitset<> *> > *POSexamples, std::vector<std::vector<boost::dynamic_bitset<> *> > *NEGexamples);
};

#endif // ONTOLOGYRULES_H
