#include "ontologyrules.h"

OntologyRules::OntologyRules(Ontology *ontoref, int ontoIndex, std::vector<Node *> *roots, std::vector< std::vector<boost::dynamic_bitset<>* > > *POSexamples, std::vector< std::vector<boost::dynamic_bitset<>* > > *NEGexamples)
{
    this->ontoref = ontoref;
    this->ontoIndex = ontoIndex;

    //example->ontology->bitset
    //std::vector< std::vector<boost::dynamic_bitset<>* > > POSexamples;
    //std::vector< std::vector<boost::dynamic_bitset<>* > > NEGexamples;
    //sthis->initBottomRules(POSexamples, NEGexamples);



    //add roots to vector
    /*for(int iroot = 0; iroot < roots->size(); ++iroot)
    {
        bottomRule newrule;
        newrule.actScore = 0;
        newrule.pathBestScore = 0;
        newrule.IDind = this->bottomRules.size();
        //covered, childs, and parents are disregarded
        this->rootsIds.push_back(newrule.IDind);
        this->bottomRules.push_back(newrule);
    }
    */
}

OntologyRules::~OntologyRules()
{

}
void OntologyRules::initBottomRules(std::vector<std::vector<boost::dynamic_bitset<> *> > *POSexamples, std::vector<std::vector<boost::dynamic_bitset<> *> > *NEGexamples)
{
 /*   //iterate over roots
    boost::dynamic_bitset<> CLOSED(this->ontoref->getOntologyParser()->getTermsNumber(), 0);
    std::vector<Node*> roots = this->ontoref->getOntologyParser()->getRoots();
    for(int iroot = 0; iroot < roots.size(); ++iroot)
    {
        std::list<Node*> OPEN;
        CLOSED[roots[iroot].idNum] = 1;
        OPEN.push_front(roots[iroot]);
        while(!OPEN.empty())
        {
            Node *top = OPEN.pop_front();
            if(CLOSED[top->idNum] == 0)
            {
                for(std::vector<edge>::iterator iedge = top->relationship.begin(); iedge != top->relationship.end(); ++iedge)
                {
                    if(iedge->end != NULL && CLOSED[iedge->end->idNum] == 0)
                    {
                        Node *newNode = iedge->end;

                        OPEN.push_back(newNode);
                    }
                }
                CLOSED[top->idNum] = 1;
            }
        }
    }


    //contains enriched terms
    this->bottomRules.resize(this->ontoref->getOntologyParser()->getTermsNumber());
    for(int irule = 0; irule < this->bottomRules.size(); ++irule)
    {
        bottomRule newRule;
        newRule.actScore = 0;
        newRule.pathBestScore = 0;
        newRule.IDind = irule;
        newRule.
    }
    */
}
