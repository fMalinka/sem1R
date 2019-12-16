#ifndef NODE_H
#define NODE_H

#include <vector>
#include <string>

struct edge;
struct termdef;

class Ontology;

//correspond to the vertice [term]
class Node
{
public:

    //variables
    std::string id; //id of the current term
    std::string name;   //the term name
    unsigned int idNum; //unique ID (in a current ontology)
    bool isRoot;    //term as a root?
    std::string def;    //definition
    std::vector<std::string> intersection_of;
    std::vector<std::string> union_of;
    std::vector<std::string> disjoint_from;
    std::vector<edge> relationship; //relationship as an edge in the graph

    bool isObsolete;    //obsolete
    std::vector<std::string> replaced_by;

    int ontoType;
    Ontology* onto_ref;

    //methods
    Node();
    ~Node();
};

//correspond to edges
struct edge
{
    //to general
    Node *start;	//pointer to the start of edge
    Node *end;	//pointer to the end of edge
    termdef *rel;
    std::string tmpend;
};

struct termdef
{
    std::string id; //id of the relationship (typedef)
    std::string name;   //the term name
    std::string def;    //definition
    bool is_cyclic;
    bool is_reflexive;
    bool is_symmetrix;
    bool is_transitive;
};
#endif // NODE_H
