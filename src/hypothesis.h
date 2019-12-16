#ifndef HYPOTHESIS_H
#define HYPOTHESIS_H

#include <vector>
#include "node.h"
#include "ontology.h"
#include <boost/unordered_map.hpp>

typedef struct {
    //boost::unordered_map<Ontology*, std::vector<Node*> > termsNode;
    std::vector<Ontology*> colOntology;
    std::vector<Ontology*> rowOntology;
    double score;
    boost::unordered_map<Ontology*,boost::unordered_map<std::string, bool> > ids;
    boost::unordered_map<Ontology*,boost::unordered_map<unsigned int, bool> > nodeIDs;  //stejne jako ids, ale key bude uint

}hypothesis;

#endif // ERROR_H
