#include "ontology.h"
Ontology::Ontology()
{;}

Ontology::Ontology(std::string name, std::string pathToOntology, Rcpp::List descriptionTerms, int ontoType)
{
    this->correct = true;
    this->ontology = NULL;

    //Check parameters
    if(name.empty())
    {
        Rcpp::Rcerr << "Fill the name of ontology!" << std::endl;
        this->correct = false;
    }
    if(pathToOntology.empty())
    {
        Rcpp::Rcerr << "Fill the path to ontology!" << std::endl;
        this->correct = false;
    }
    if(descriptionTerms.size() == 0)
    {
        Rcpp::Rcerr << "Set description of terms!" << std::endl;
        this->correct = false;
    }

    if(this->correct)
    {
        this->name = name;
        this->pathToOntology = pathToOntology;
        this->descriptionTerms = convertRList2Vector(descriptionTerms);
        this->ontoType = ontoType;

        Rcpp::Rcout << "[" << name  << " ONTOLOGY]" << std::endl;
        OntologyParser *onto = new OntologyParser(pathToOntology.c_str(), this->ontoType, this);
        int state = onto->loadFile();
        if(state == OK)
        {
            this->ontology = onto;
            Rcpp::Rcout << "Number of Terms: " << onto->getTermsNumber() << " (NON-obsolete)" << std::endl;
            onto->printRoots();
            //ontoRow.printTermIds();
            //onto->printRoots();
        }
        else
        {
            //if(state == FILE_NOT_FOUND)
            //    Rcpp::Rcerr << "Ontology file not found!" << std::endl;
            this->correct = false;
        }
    }

}

Ontology::~Ontology()
{
    for(int i = 0; i < this->SemanticPatterns.size(); ++i)
    {
        delete this->SemanticPatterns[i];
    }
}

int Ontology::setCoverageVector(hypothesis *searched, arma::vec *bitset, int *arrayIndex)
{
    //initialize - find nodes for corresponding terms
    std::list<Node*> nodePool;
    std::vector<std::string> terms = this->descriptionTerms[*arrayIndex];
    for(std::vector<std::string>::iterator iterm = terms.begin(); iterm != terms.end(); ++iterm)
    {
        Node *termNode = this->getOntologyParser()->getNodesBottomUP(&(*iterm));
        if(termNode != NULL)
             nodePool.push_front(termNode);
    }

    //denote corresponding element
    if(findTestAncestors2(&nodePool, searched) == OK)
        (*bitset)(*arrayIndex) = 1;

    return OK;
}

void Ontology::printSemanticPattern(boost::dynamic_bitset<> *bitset)
{
    for(int ibit = 0; ibit < bitset->size(); ++ibit)
    {
        if((*bitset)[ibit] == 1)
        {
            Rcpp::Rcout << this->getOntologyParser()->getNodesByPosition(ibit)->name << "|";
        }
    }
    Rcpp::Rcout << std::endl;
}

int Ontology::isCovered(Node *searched, arma::vec *bitset, int *arrayIndex)
{
    //initialize - find nodes for corresponding terms
    std::list<Node*> nodePool;
    std::vector<std::string> terms = this->descriptionTerms[*arrayIndex];
    for(std::vector<std::string>::iterator iterm = terms.begin(); iterm != terms.end(); ++iterm)
    {
        Node *termNode = this->getOntologyParser()->getNodesBottomUP(&(*iterm));
        if(termNode != NULL)
             nodePool.push_front(termNode);
    }

    //denote corresponding element
    if(findTestAncestors(&nodePool, searched) == OK)
        (*bitset)(*arrayIndex) = 1;

    return OK;
}

int Ontology::findTestAncestors2(std::list<Node*> *nodePool, hypothesis *searched)
{
    if(nodePool->empty())
        return TERM_NOT_FOUND;
    //we use list -> LIFO -> DFS
    //1.push the first element
    Node *front = nodePool->front();

    //2.find in ontology, test for equivalence
    if(searched->ids[this].count(front->id)) //TODO -> rychlejsi bude neporovnavat string!
        return OK;

    boost::unordered_map<std::string, bool> CLOSED;
    while(!nodePool->empty())
    {
        //1.push the first element
        front = nodePool->front();
        nodePool->pop_front();
        //3. if node is not equivalent, generate new ancestors and add them at the beggining
        for(std::vector<edge>::iterator iedge = front->relationship.begin(); iedge != front->relationship.end(); ++iedge)
        {
            if(iedge->end != NULL)
            {
                //b) nasel jsi, nastav covered = true
                if(searched->ids[this].count(iedge->end->id))
                    return OK;
                if(CLOSED.find(iedge->end->id) == CLOSED.end())
                    nodePool->push_front(iedge->end);
            }
        }
    }
    return TERM_NOT_FOUND;
}


int Ontology::findTestAncestors(std::list<Node*> *nodePool, Node *searched)
{
    if(nodePool->empty())
        return TERM_NOT_FOUND;
    //we use list -> LIFO -> DFS
    //1.push the first element
    Node *front = nodePool->front();

    //2.find in ontology, test for equivalence
    if(front->id == searched->id) //TODO -> rychlejsi bude neporovnavat string!
        return OK;

    boost::unordered_map<std::string, bool> CLOSED;
    while(!nodePool->empty())
    {
        //1.push the first element
        front = nodePool->front();
        nodePool->pop_front();
        //3. if node is not equivalent, generate new ancestors and add them at the beggining
        for(std::vector<edge>::iterator iedge = front->relationship.begin(); iedge != front->relationship.end(); ++iedge)
        {
            if(iedge->end != NULL)
            {
                //b) nasel jsi, nastav covered = true
                if(iedge->end->id == searched->id)
                    return OK;
                if(CLOSED.find(iedge->end->id) == CLOSED.end())
                    nodePool->push_front(iedge->end);
            }
        }
    }
    return TERM_NOT_FOUND;
}

int Ontology::getLevelOfSpecialization(Node *specific, Node *general)
{
    int level = 0;
    std::list<Node*> nodePool;
    nodePool.push_back(specific);
    //we use list -> FIFO -> BFS
    //1.push the first element
    Node *front = nodePool.front();

    //2.find in ontology, test for equivalence
    if(front->id == general->id) //TODO -> rychlejsi bude neporovnavat string!
        return level;

    boost::unordered_map<std::string, bool> CLOSED;
    while(!nodePool.empty())
    {
        ++level;
        //1.push the first element
        front = nodePool.front();
        nodePool.pop_front();
        //3. if node is not equivalent, generate new ancestors and add them at the beggining
        for(std::vector<edge>::iterator iedge = front->relationship.begin(); iedge != front->relationship.end(); ++iedge)
        {
            if(iedge->end != NULL)
            {
                //b) nasel jsi, nastav covered = true
                if(iedge->end->id == general->id)
                    return level;
                if(CLOSED.find(iedge->end->id) == CLOSED.end())
                    nodePool.push_back(iedge->end);
            }
        }
    }
    return TERM_NOT_FOUND;
}


int Ontology::findLGNancestors(Node *searched, boost::dynamic_bitset<> *bitset, boost::unordered_map<Node*, int> *mapNode2Int)
{
    std::list<Node*> nodePool;
    nodePool.push_front(searched);

    boost::unordered_map<std::string, bool> CLOSED;
    while(!nodePool.empty())
    {
        Node *front = nodePool.front();
        nodePool.pop_front();

        //record
       bitset->set((*mapNode2Int)[front], 1);//(*bitset)[(*mapNode2Int)[front]] = 1;

        //3. if node is not equivalent, generate new ancestors and add them at the beggining
        for(std::vector<edge>::iterator iedge = front->relationship.begin(); iedge != front->relationship.end(); ++iedge)
        {
            if(iedge->end != NULL)
            {
                if(CLOSED.find(iedge->end->id) == CLOSED.end())
                    nodePool.push_front(iedge->end);
            }
        }
    }
    return OK;
}

boost::dynamic_bitset<> Ontology::getTermBitAncestors(Node *searched)
{
    std::list<Node*> nodePool;
    nodePool.push_front(searched);

    boost::dynamic_bitset<> bitset(this->getOntologyParser()->getTermsNumber(), 0);
    //boost::dynamic_bitset<> CLOSED(this->getOntologyParser()->getTermsNumber(), 0);

    //boost::unordered_map<std::string, bool> CLOSED;
    while(!nodePool.empty())
    {
        Node *front = nodePool.front();
        nodePool.pop_front();

        //record 1
       bitset.set(front->idNum, 1);

        //3. if node is not equivalent, generate new ancestors and add them at the beggining
        for(std::vector<edge>::iterator iedge = front->relationship.begin(); iedge != front->relationship.end(); ++iedge)
        {
            if(iedge->end != NULL)
            {
                //if(CLOSED.find(iedge->end->id) == CLOSED.end())
                //    nodePool.push_front(iedge->end);
                if(bitset[iedge->end->idNum] == 0)
                    nodePool.push_back(iedge->end);
            }
        }
    }
    return bitset;
}

std::vector<Node*> Ontology::findLGNroots()
{
    //vytvor hasmapu, kde kazdemu termu v ontologii je prirazeno cislo reprezentujici pozici v bitset ci matici
    //bitset i = 1 pro termu ktery jsou rodicem
    //vsechny bitsety z-andovat, coz reprezentuje spolecnou podmnozinu termu pro cely dataset -> zkusit udelat pomoci matice

    //pro nalezeni LGN je postup:
        //prochazej ontologii top-down, pokud je na prislusnym miste 1, pokracuj dal ... pokud uz neni, vezmi posledni


    std::vector<Node*> LGN(this->getOntologyParser()->getRoots()->size());
    std::vector<std::vector<std::string> > separatedSets(this->getOntologyParser()->getRoots()->size());
    //separates description terms due to roots
    #pragma omp parallel for
    for(int iroot = 0; iroot < this->getOntologyParser()->getRoots()->size(); ++iroot)
    {
        for(std::vector<std::vector<std::string> >::iterator idesc = this->descriptionTerms.begin(); idesc != this->descriptionTerms.end(); ++idesc)
        {
            for(std::vector<std::string>::iterator iidesc = idesc->begin(); iidesc != idesc->end(); ++iidesc)
            {
                std::list<Node*> nodePool;
                nodePool.push_front((*this->getOntologyParser()->getBottomUP())[*iidesc]);
                if(findTestAncestors(&nodePool, (*this->getOntologyParser()->getRoots())[iroot]) == OK)
                    separatedSets[iroot].push_back(*iidesc);
            }
        }
    }

    //make a hashes
    std::vector<Node*> mapInt2Node;
    boost::unordered_map<Node*, int> mapNode2Int;
    int position = 0;
    boost::unordered_map<std::string, Node*> *ontologyMap = this->getOntologyParser()->getBottomUP();
    for(boost::unordered_map<std::string, Node*>::iterator imap = ontologyMap->begin(); imap != ontologyMap->end(); ++imap)
    {
        mapNode2Int[imap->second] = position;
        mapInt2Node.push_back(imap->second);
        ++position;
    }

    //separates description terms due to roots
    #pragma omp parallel for
    for(int iroot = 0; iroot < this->getOntologyParser()->getRoots()->size(); ++iroot)
    {
        boost::dynamic_bitset<> bitsetLGN(this->getOntologyParser()->getTermsNumber(), 1);

        bitsetLGN.set();
        //iterate over all nodes
        for(int i = 0; i < separatedSets[iroot].size(); ++i)
        {
            boost::dynamic_bitset<> bit(this->getOntologyParser()->getTermsNumber(), 0);
            this->findLGNancestors((*this->getOntologyParser()->getBottomUP())[separatedSets[iroot][i]], &bit, &mapNode2Int);
            //AND operation, we try to find an intersect
            bitsetLGN &= bit;
        }

        //get name of nodes
        int val_max = 0;
        Node *node_max = NULL;
        for(int i = 0; i < bitsetLGN.size(); ++i)
        {
            if(bitsetLGN.test(i) == 1)//if(bitsetLGN[i] == 1)
            {
                Node *specific = (*this->getOntologyParser()->getBottomUP())[mapInt2Node[i]->id];
                Node *root = (*this->getOntologyParser()->getBottomUP())[(*this->getOntologyParser()->getRoots())[iroot]->id];
                int level = this->getLevelOfSpecialization(specific, root);

                //choose a node with maximum value (more specific)
                if(level >= val_max)
                {
                    val_max = level;
                    node_max = (*this->getOntologyParser()->getTopDown())[mapInt2Node[i]->id];
                }
            }
        }
        LGN[iroot] = node_max;
    }

return LGN;
}

std::vector<std::vector<std::string> > Ontology::convertRList2Vector(Rcpp::List desc)
{
    int n = desc.size();
    std::vector<std::vector<std::string> > newDesc(n);
    for(int i = 0; i < n; ++i)
    {
        Rcpp::StringVector stringVector = desc[i];
        std::vector<std::string> tmp_desc;
        for(int ii = 0; ii < stringVector.size(); ++ii)
        {
            tmp_desc.push_back(Rcpp::as<std::string>(stringVector[ii]));
        }
        newDesc[i] = tmp_desc;
    }
    return newDesc;
}

void Ontology::precomputeSemanticPatterns(int vectorsize)
{
    this->SemanticPatterns.resize(vectorsize);
    for(int ivec = 0; ivec < vectorsize; ++ivec)
    {
        boost::dynamic_bitset<> *pattern = new boost::dynamic_bitset<>(this->getOntologyParser()->getTermsNumber(),0);

        //terms from descriptions
        std::vector<std::string> terms = this->descriptionTerms[ivec];

        for(std::vector<std::string>::iterator iterm = terms.begin(); iterm != terms.end(); ++iterm)
        {
            Node *termNode = this->getOntologyParser()->getNodesBottomUP(&(*iterm));
            if(termNode != NULL)
            {
                *pattern |= this->getTermBitAncestors(termNode);
            }
        }
        //std::cout << "pattern: " << ivec << " count:" << pattern->count() << std::endl;
        this->SemanticPatterns[ivec] = pattern;
    }
/*
    //matrix initialization
    this->precomputedDistanceMatrix.resize(vectorsize);
    for(int i = 0; i < vectorsize; ++i)
    {
        this->precomputedDistanceMatrix[i].resize(vectorsize);
    }

    //precompute score as a distance matrix
    for(int r = 0; r < vectorsize; ++r)
    {
        for(int c = r+1; c < vectorsize; ++c)
        {
            int common = ((*SemanticPatterns[r]) | (*SemanticPatterns[c])).count();
            int intersect = ((*SemanticPatterns[r]) & (*SemanticPatterns[c])).count();

            double score = 0;
            if(common != 0)
            {
                score = intersect/(double)common;
            }
            //symetric relationship
            this->precomputedDistanceMatrix[r][c] = score;
            this->precomputedDistanceMatrix[c][r] = score;
        }
    }
    */
}

void Ontology::precomputeSemanticPatternsTest(int vectorsize)
{
    this->SemanticPatternsTest.resize(vectorsize);
    for(int ivec = 0; ivec < vectorsize; ++ivec)
    {
        boost::dynamic_bitset<> *pattern = new boost::dynamic_bitset<>(this->getOntologyParser()->getTermsNumber(),0);

        //terms from descriptions
        std::vector<std::string> terms = this->descriptionTermsTest[ivec];

        for(std::vector<std::string>::iterator iterm = terms.begin(); iterm != terms.end(); ++iterm)
        {
            Node *termNode = this->getOntologyParser()->getNodesBottomUP(&(*iterm));
            if(termNode != NULL)
            {
                *pattern |= this->getTermBitAncestors(termNode);
            }
        }
        //std::cout << "pattern: " << ivec << " count:" << pattern->count() << std::endl;
        this->SemanticPatternsTest[ivec] = pattern;
    }
}
