#ifndef ONTOLOGY_PARSER_H
#define ONTOLOGY_PARSER_H

#include <iostream>
#include <fstream>
#include <string>
#include <boost/regex.hpp>
#include <boost/unordered_map.hpp>
#include <vector>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <RcppArmadillo.h>

#include "node.h"
#include "error.h"

class Ontology;


class OntologyParser {
public:
    OntologyParser(const char *path, int ontoType, Ontology* onto_ref): path(path), ontoType(ontoType), onto_ref(onto_ref) {}
    int loadFile();
    void printTermIds();
    void printMoreGeneral(std::string id);
    void printMoreSpecific(std::string id);
    void printRoots();
    int getTermsNumber();
    std::vector<Node*>* getRoots(){return &(this->roots);}
    Node* getNodesBottomUP(std::string *name) {return this->mapNodesBottomUP[*name];}
    Node* getNodeTopDown(std::string *name) {return this->mapNodesTopDown[*name];}
    boost::unordered_map<std::string, Node*>* getBottomUP() {return &(this->mapNodesBottomUP);}
    boost::unordered_map<std::string, Node*>* getTopDown() {return &(this->mapNodesTopDown);}
    Node* getNodesByPosition(int pos){return this->Nodes_topDown[pos];}
    Node* getNodesByPosition_bottomUP(int pos){return this->Nodes_bottomUP[pos];}

private:
    const char *path;
    std::vector<Node*> Nodes;
    std::vector<Node*> Nodes_bottomUP;
    std::vector<Node*> Nodes_topDown;
    boost::unordered_map<std::string, Node*> mapNodesBottomUP;
    boost::unordered_map<std::string, Node*> mapNodesTopDown;
    std::vector<Node*> roots;
    int ontoType;

    Ontology* onto_ref;

    void parseData(std::vector<std::string> *dataTerm , std::vector<std::string> *dataTypedef);
    void loadItem(std::vector<std::string> *data, std::ifstream *myfile);
    struct termdef createTermDef(std::string *text);
    Node* createTerm(std::string *text, std::vector<termdef> *typedefs, unsigned int *index);
    std::vector<Node*> makeTopDownNodes(std::vector<Node*> *nodes);
    bool checkBottomTop(std::vector<std::string> *keys);
    std::vector<Node*> findRoots(std::vector<std::string> *allKeys);

};

#endif


