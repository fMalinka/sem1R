#include "ontologyParser.h"
#include <omp.h>


void OntologyParser::loadItem(std::vector<std::string> *data, std::ifstream *myfile)
{
    std::string line = "";
    std::string text = "#";
    while(getline(*myfile, line))
    {
        text += line;
        text += "#";
        if(line.empty())    //separate by a new line (blank line)
        {
            data->push_back(text);
            return;
        }
    }
}

struct termdef OntologyParser::createTermDef(std::string *text)
{
    boost::match_results<std::string::const_iterator> res;
    struct termdef myterm;  //my new struct
    if(boost::regex_search(*text,res, boost::regex("#id:\\s*([^\\s]+)#", boost::regex::icase)))
    {
        myterm.id = res[1].str();
    }
    if(boost::regex_search(*text,res, boost::regex("#name:\\s*([^\\#]+)", boost::regex::icase)))
    {
        myterm.name = res[1].str();
    }
    if(boost::regex_search(*text,res, boost::regex("#id:\\s*([^\\#]+)", boost::regex::icase)))
    {
        myterm.def = res[1].str();
    }
    if(boost::regex_search(*text,res, boost::regex("#is_cyclic:\\s*true", boost::regex::icase)))
    {
        myterm.is_cyclic = true;
    }
    else
        myterm.is_cyclic = false;
    if(boost::regex_search(*text,res, boost::regex("#is_reflexive:\\s*true", boost::regex::icase)))
    {
        myterm.is_reflexive = true;
    }
    else
        myterm.is_reflexive = false;
    if(boost::regex_search(*text,res, boost::regex("#is_symmetrix:\\s*true", boost::regex::icase)))
    {
        myterm.is_symmetrix = true;
    }
    else
        myterm.is_symmetrix = false;
    if(boost::regex_search(*text,res, boost::regex("#is_transitive:\\s*true", boost::regex::icase)))
    {
        myterm.is_transitive = true;
    }
    else
        myterm.is_transitive = false;

    return myterm;
}

Node* OntologyParser::createTerm(std::string *text, std::vector<termdef> *typedefs, unsigned int *index)
{
    boost::match_results<std::string::const_iterator> res;
    //Node *mynode = new Node();  //my new node
    //mynode->isObsolete = false;
    if(boost::regex_search(*text,res, boost::regex("#is_obsolete:\\s*true", boost::regex::icase)))
    {
        return NULL;
    }
    else
    {
        Node *mynode = new Node();  //my new node
        mynode->isObsolete = false;

        if(boost::regex_search(*text,res, boost::regex("#id:\\s*([^\\s]+)#", boost::regex::icase)))
        {
            mynode->id = res[1].str();
            //std::cout << res[1].str() << std::endl;
        }
        if(boost::regex_search(*text,res, boost::regex("#name:\\s*([^\\#]+)", boost::regex::icase)))
        {
            mynode->name = res[1].str();
        }
        if(boost::regex_search(*text,res, boost::regex("#def:\\s*([^\\#]+)", boost::regex::icase)))
        {
            mynode->def = res[1].str();
        }

        //type of ontology row/column
        mynode->ontoType = this->ontoType;
        mynode->onto_ref = this->onto_ref;
        //mynode->idNum = *index;

        for(std::vector<termdef>::iterator idef = typedefs->begin(); idef != typedefs->end(); ++idef)
        {
    /*
      old version
            boost::match_results<std::string::const_iterator> what;
            if(boost::regex_search(*text,what, boost::regex(std::string((*idef).id + ":?\\s*([^\\s]+)"), boost::regex::icase)))
            {
                struct edge e = {.start = mynode, .end = NULL, .rel = &(*idef), .tmpend = what[1].str()};  //ukazatel na struct????
                mynode->relationship.push_back(e);
            }
    */
            //std::cout << (*idef).id << std::endl;
            boost::regex expr = boost::regex(std::string("#(relationship:\\s*)?" + (*idef).id + ":?\\s*([^\\s]+)"), boost::regex::icase);
            boost::regex_token_iterator<std::string::iterator> iter(text->begin(), text->end(), expr, 2);
            boost::regex_token_iterator<std::string::iterator> end;
            for(;iter != end; ++iter )
            {
                struct edge e = {.start = mynode, .end = NULL, .rel = &(*idef), .tmpend = *iter};  //ukazatel na struct????
                mynode->relationship.push_back(e);
                //std::cout << (*idef).id <<  " "  << *iter << std::endl;
            }

        }
        return mynode;
    }
    return NULL;
}


int OntologyParser::loadFile()
{
        std::ifstream myfile(this->path);
	std::vector<std::string> dataTerm;
	std::vector<std::string> dataTypedef;
	if(myfile.is_open())
	{
            std::string line;
            boost::smatch what;
            Rcpp::Rcout << "Loading \"" << this->path << "\" ...";
            while(getline(myfile, line))
            {
                if(boost::regex_match(line, what, boost::regex("\\s*\\[Term\\]", boost::regex::icase)))
                {
                loadItem(&dataTerm, &myfile);
                }
                else if(boost::regex_match(line, what, boost::regex("\\s*\\[Typedef\\]", boost::regex::icase)))
                {
                loadItem(&dataTypedef, &myfile);
                }
            }

            parseData(&dataTerm , &dataTypedef);
            Rcpp::Rcout << " DONE" << std::endl;
	}
	else
	{
                Rcpp::Rcerr << "Cannot open: \"" << this->path << "\"" << std::endl;
                return FILE_NOT_FOUND;
	}
        myfile.close();
    return OK;
}

void OntologyParser::parseData(std::vector<std::string> *dataTerm , std::vector<std::string> *dataTypedef)
{
    //parse Typedef
    int sizeTypedef = dataTypedef->size();
    std::vector<termdef> typedefs(sizeTypedef);
    #pragma omp parallel for
    for(int i = 0; i < sizeTypedef; ++i)
    {
        typedefs[i] = createTermDef(&(*dataTypedef)[i]);
    }
    struct termdef is_a = {.id = "is_a", .name="is a", .def = "", .is_cyclic = true, .is_reflexive = true, .is_symmetrix = true, .is_transitive = true};
    typedefs.push_back(is_a);
    sizeTypedef = dataTypedef->size();

    //parse Node
    int sizeNode = dataTerm->size();
    std::vector<Node*> typeNodes(sizeNode);

    #pragma omp parallel for
    for(unsigned int i = 0; i < sizeNode; ++i)
    {
        Node *newTerm = createTerm(&(*dataTerm)[i], &typedefs, &i);
        if(newTerm != NULL)
            typeNodes[i] = newTerm;
    }
    //update the size of vector --- try to do with OPENMP
    std::vector<Node*> tmpNodes = typeNodes;
    unsigned int inonObsolete = 0;
    typeNodes.clear();
    for(int i = 0; i < tmpNodes.size(); ++i)
    {
        if(tmpNodes[i] != NULL)
        {
            tmpNodes[i]->idNum = inonObsolete;
            typeNodes.push_back(tmpNodes[i]);
            inonObsolete++;
        }
    }
    sizeNode = typeNodes.size();
    std::vector<std::string> keys;
    //make a hash
    for(int i = 0; i < sizeNode; ++i)
    {
        this->mapNodesBottomUP[typeNodes.at(i)->id] = typeNodes.at(i);
        keys.push_back(typeNodes.at(i)->id);
    }
    #pragma omp parallel for
    for(int i = 0; i < sizeNode; ++i)
    {
        for(std::vector<edge>::iterator iedge = typeNodes.at(i)->relationship.begin(); iedge != typeNodes.at(i)->relationship.end(); ++iedge)
        {
            if(this->mapNodesBottomUP.find((*iedge).tmpend) != this->mapNodesBottomUP.end())
                (*iedge).end = this->mapNodesBottomUP[(*iedge).tmpend];
            else
                Rcpp::Rcerr << std::endl << "Relationship '" << (*iedge).tmpend << "' in term '" << typeNodes.at(i)->id << "' is currupted! [IGNORED]";
            //to specific
        }
    }
    this->Nodes_bottomUP = typeNodes;
    this->Nodes = typeNodes;
    this->Nodes_topDown = this->makeTopDownNodes(&typeNodes);
    
  /*
    for(int i = 0; i < this->Nodes.size(); ++i)
    {
		Node *mynode = this->Nodes[i];
		std::cout << "id:#" << mynode->id << "#" << std::endl;
		std::cout << "name:#" << mynode->name << "#" << std::endl;;
		std::cout << "def:#" << mynode->def << "#" << std::endl;;
		std::vector<edge> relationship = mynode->relationship;
		for(int ii = 0; ii < relationship.size(); ++ii)
		{
			if(relationship[ii].end != NULL)
			{
				std::cout << "is_a:#" << relationship[ii].end->id << "#" << std::endl;
			}
			
		}
		std::cout << std::endl << std::endl;
	}
		//std::vector<Node*>
*/
    //find roots
    this->roots = findRoots(&keys);


    //Equality test
    //bool test = checkBottomTop(&keys);
    //std::cout << "Bottom-UP and TopDown are equal: " << test << std::endl;
}

void OntologyParser::printTermIds()
{
    for(std::vector<Node*>::iterator it = this->Nodes.begin(); it != this->Nodes.end(); ++it)
    {
        if((*it) != NULL)
            Rcpp::Rcout << (*it)->id << std::endl;
    }
    Rcpp::Rcout << "Number of terms: " << this->Nodes.size() << std::endl;
}

int OntologyParser::getTermsNumber()
{
    return this->Nodes.size();
}

std::vector<Node*> OntologyParser::makeTopDownNodes(std::vector<Node*> *nodes)
{
    //dont use OPENMP
    std::vector<Node*> noNullNodes;
    std::vector<Node*> nodesTopDown;
    boost::unordered_map<std::string, Node*> mapNodes;
    for(std::vector<Node*>::iterator it = nodes->begin(); it != nodes->end(); ++it)
    {
        //omit null nodes (ie. obsolete) and construct new nodes
        if((*it) != NULL)
        {
            Node *newNode = new Node(*(*it));   //implicit copy constructor
            Node *newNodeTopDown = new Node(*(*it));   //implicit copy constructor
            newNodeTopDown->relationship.clear();
            noNullNodes.push_back(newNode);
            mapNodes[newNodeTopDown->id] = newNodeTopDown;
            nodesTopDown.push_back(newNodeTopDown);
        }
    }

    for(int i = 0; i < noNullNodes.size(); ++i)//std::vector<Node*>::iterator it = noNullNodes.begin(); it != noNullNodes.end(); ++it)
    {

            //go to the parents
            for(std::vector<edge>::iterator itEdge = noNullNodes[i]->relationship.begin(); itEdge != noNullNodes[i]->relationship.end(); ++itEdge)
            {

                if((*itEdge).end != NULL)
                {
                    Node *parent = (*itEdge).end;
                    edge newEdge = (*itEdge);
                    newEdge.start = mapNodes[parent->id];
                    newEdge.end = nodesTopDown[i];
                    mapNodes[parent->id]->relationship.push_back(newEdge);
                }

            }

    }

    //clean up
    for(int i = 0; i < noNullNodes.size(); ++i)
    {
            delete noNullNodes[i];
    }
    this->mapNodesTopDown = mapNodes;
    return nodesTopDown;
}

bool OntologyParser::checkBottomTop(std::vector<std::string> *keys)
{
    for(std::vector<std::string>::iterator it = keys->begin(); it != keys->end(); ++it)
    {
        Node *tmpNode= this->mapNodesBottomUP[(*it)];

        if(tmpNode == NULL)
        {
           Rcpp::Rcout << "NULL: " << *it << std::endl;
            continue;
        }
        //std::cout << "check: " << tmpNode->id << std::endl;
        for(std::vector<edge>::iterator itEdge = tmpNode->relationship.begin(); itEdge != tmpNode->relationship.end(); ++itEdge)
        {
            Node *top = itEdge->end;
            if(top == NULL)
                continue;
            Node *topTopDown = this->mapNodesTopDown[top->id];
            bool same = false;
            for(std::vector<edge>::iterator itEdgeTop = topTopDown->relationship.begin(); itEdgeTop != topTopDown->relationship.end(); ++itEdgeTop)
            {
                if(itEdgeTop->end->id == tmpNode->id)
                {
                    same = true;
                }
            }
            if(same == false)
                return false;
        }
    }
    return true;
}

void OntologyParser::printMoreGeneral(std::string id)
{
    Node *actNode = this->mapNodesBottomUP[id];    
    if(actNode == NULL)
        Rcpp::Rcout << id  << " NOT FOUND" << std::endl;
    else
    {
        Rcpp::Rcout << "More general from (" << actNode->relationship.size() << ") "  << id << " is:";
        for(std::vector<edge>::iterator it = actNode->relationship.begin(); it != actNode->relationship.end(); ++it)
        {
            if(it->end != NULL)
                Rcpp::Rcout << " " << it->end->id;
        }
    }
    Rcpp::Rcout << std::endl;
}

void OntologyParser::printMoreSpecific(std::string id)
{
    Node *actNode = this->mapNodesTopDown[id];
    if(actNode == NULL)
        Rcpp::Rcout << id  << " NOT FOUND" << std::endl;
    else
    {
        Rcpp::Rcout << "More specific from (" << actNode->relationship.size() << ") "  << id << " is:";
        for(std::vector<edge>::iterator it = actNode->relationship.begin(); it != actNode->relationship.end(); ++it)
        {
            if(it->end != NULL)
                Rcpp::Rcout << " " << it->end->id;
        }
    }
    Rcpp::Rcout << std::endl;
}

std::vector<Node*> OntologyParser::findRoots(std::vector<std::string> *allKeys)
{
    std::vector<Node*> roots;
    //TODO: uvazuj o OPENMP -> v kazdym vytvorim subpole, ktery pak sloucim do jednoho
    for(std::vector<std::string>::iterator iKey = allKeys->begin(); iKey != allKeys->end(); ++iKey)
    {
        Node *actNode = this->mapNodesBottomUP[*iKey];
        if(actNode->relationship.empty())
        {
            roots.push_back(actNode);
        }
    }
    return roots;
}

void OntologyParser::printRoots()
{
    Rcpp::Rcout << "Roots:";
    for(std::vector<Node*>::iterator it = this->roots.begin(); it != this->roots.end(); ++it)
    {
        Rcpp::Rcout << " " << (*it)->id;
    }
    Rcpp::Rcout << std::endl;
}
