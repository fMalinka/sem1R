#include "conjunction_learning.h"

#define EPS_DOUBLE 0.01

//Global parametrs with DEFAULTs values
string TRAIN;
bool TRAIN_FLAG = 0;
//string TEST;
bool TEST_FLAG = 0;
int QUEUE_LIMIT = 100;
bool QUEUE_UNLIMITED = false;
int RULES_LIMIT = 10;
int ITERATE_LIMIT = 100000;
int VERBOSE = 0;
int DEBUG = 0;
int ROC_FLAG = 0;
//string ROC;

//HASH FOR BITSET https://stackoverflow.com/questions/45914271/generate-hash-for-boostdynamic-bitset-and-convert-hash-back-to-boostdynamic
//https://stackoverflow.com/questions/3896357/unordered-hash-map-from-bitset-to-bitset-on-boost
//#define BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS
//namespace boost {
//    template <typename B, typename A>
//    std::size_t hash_value(const boost::dynamic_bitset<B, A>& bs) {
//        return boost::hash_value(bs.m_bits);
//    }
//}
/*
#define BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS
namespace boost {
    template <typename B, typename A>
    inline std::size_t hash_value(const boost::dynamic_bitset<B, A>& bs) {
        std::vector<B, A> v;
        boost::to_block_range(bs, std::back_inserter(v));
        return boost::hash_value(v);
    }
}
*/

BeamTopDown::BeamTopDown(arma::mat *data, std::vector<Ontology *> *refOntologies, std::vector< std::vector<boost::dynamic_bitset<>* > > *POSexamples, std::vector< std::vector<boost::dynamic_bitset<>* > > *NEGexamples, std::vector<boost::dynamic_bitset<> > *enrichItems)
{
    this->data = data;
    this->refOntologies = refOntologies;
    this->POSexamples = POSexamples;
    this->NEGexamples = NEGexamples;
    this->enrichItems = enrichItems;
}
BeamTopDown::BeamTopDown(std::vector<bottomFeature> *bottomFeatures, boost::dynamic_bitset<> *myclass, std::string objective, int verbose)
{
    this->bottomFeatures = bottomFeatures;
    this->myclass = myclass;
    this->verbose = verbose;
     
    if(objective == "f1")
    {
		evaluationFunction = &BeamTopDown::evaluateF1Score;
		evaluationFunctionPotential = &BeamTopDown::evaluateF1ScorePotencialBoth;
	}
	else if(objective == "acc")
	{
		evaluationFunction = &BeamTopDown::evaluateACCScore;
		evaluationFunctionPotential = &BeamTopDown::evaluateACCScorePotencialBoth;
	}
	else if(objective == "auc")
	{
		evaluationFunction = &BeamTopDown::evaluateAUCcore;
		evaluationFunctionPotential = &BeamTopDown::evaluateAUCScorePotencialBoth;
	}
	else
	{
                Rcpp::Rcout << "Objective function: " << objective << " has not been found! F1 is set as an imlicit one." << std::endl;
		evaluationFunction = &BeamTopDown::evaluateF1Score;
		evaluationFunctionPotential = &BeamTopDown::evaluateF1ScorePotencialBoth;
	}
}

newComplex BeamTopDown::runExhaustive(int filterTH, int ruleDepth, double signTH)
{
    //first step - initialize

    clock_t begin = std::clock();

    int iterate = 0;
    Rcpp::Rcout << "runExhaustive: " << iterate << std::endl;
    newComplex bestRule;
    bestRule.score = 0;
    bool atleastOneGenerated = false;

    //newComplex bestRuleEntropy;
    //bestRuleEntropy.score = INF;

    //newComplex bestRuleACC;
    //bestRuleACC.score = 0;

    //std::list<newComplex> newCandidates;
    std::priority_queue<newComplex> newCandidates;
    for(int ifeature = 0; ifeature < this->bottomFeatures->size(); ++ifeature)
    {
        newComplex newItem;
        //newItem.score = evaluateF1Score(&((*bottomFeatures)[ifeature].exampleCovered), this->myclass);
        //newItem.score = evaluateACCScore(&((*bottomFeatures)[ifeature].exampleCovered), this->myclass);
        newItem.score = (this->*evaluationFunction)(&((*bottomFeatures)[ifeature].exampleCovered), this->myclass);
        newItem.rules.push_back(&((*bottomFeatures)[ifeature]));
        newCandidates.push(newItem);

//       if(bestRule.score <= newItem.score)
//            bestRule = newItem;


		double sig = evaluateSignificance(&newItem, myclass);               
        newItem.sig = sig;
        
        newCandidates.push(newItem);
        //std::cout << "ssig: " << sig << std::endl;
     //   if(sig > signTH)
     //   {
			atleastOneGenerated = true;
            if(bestRule.score <= newItem.score)
                bestRule = newItem;
     //   }

        //std::cout << "IDint: " << newItem.rules[0]->IDind <<  " ";

        //double entropy = evaluateEntropyScore(&((*bottomFeatures)[ifeature].exampleCovered), this->myclass);
        //if(bestRuleEntropy.score > entropy)
        //{
        //    bestRuleEntropy = newItem;
        //    bestRuleEntropy.score = entropy;
        //}

        //double acc = evaluateACCScore(&((*bottomFeatures)[ifeature].exampleCovered), this->myclass);
        //if(bestRuleACC.score < acc)
        //{
        //    bestRuleACC = newItem;
        //    bestRuleACC.score = acc;
        //}

        ++iterate;
    }
    Rcpp::Rcout << std::endl;
    //Rcpp::Rcout << "runExhaustive: " << iterate << " best score: " << bestRule.score << " rule: " << getPrintableFeature(&bestRule) << std::endl;

    Rcpp::Rcout << "RUN 1 - size: " << iterate << " best score: " << bestRule.score << std::endl;
    Rcpp::Rcout << "Rule: " << getPrintableFeature(&bestRule) << std::endl << std::endl;
    clock_t end = std::clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    Rcpp::Rcout << "elapsed time: " << elapsed_secs << std::endl << std::endl;
    filterComplexes(&newCandidates, filterTH);
    //if(!atleastOneGenerated)
    //{
	//	return bestRule;
	//}
    //iterate
    for(int icycle = 0; icycle < ruleDepth; ++icycle)
    {
        iterate = 0;
        atleastOneGenerated = false;
        std::priority_queue<newComplex> newTmp;        
        boost::unordered_map<boost::dynamic_bitset<>, bool> duplicitiesBITSET;
        //bestRule.score = 0;
        int iii = 0;
        while(!newCandidates.empty())
        {
            newComplex toExpand = newCandidates.top();
            newCandidates.pop();
            duplicitiesBITSET[getRuleBitset(&toExpand)] = true;
            for(int ipos = 0; ipos < this->bottomFeatures->size(); ++ipos)
            {
                bottomFeature *mycandFeature = &((*bottomFeatures)[ipos]);
                newComplex newItem;
                newItem.rules = toExpand.rules;
                newItem.rules.push_back(mycandFeature);
                if(duplicitiesBITSET.find(getRuleBitset(&newItem)) == duplicitiesBITSET.end())
                {
                    boost::dynamic_bitset<> cover(this->myclass->size());
                    cover.set();
                    bool hasColOntology = false;
                    for(int irule = 0; irule < newItem.rules.size(); ++irule)
                    {
                                                Rcpp::Rcout << " b: " << newItem.rules[irule]->exampleCovered.count();
                        cover = cover & newItem.rules[irule]->exampleCovered;
                        if(newItem.rules[irule]->noderef->ontoType != ROW_ONTOLOGY)
							hasColOntology = true;
                    }
                    Rcpp::Rcout << std::endl;
                    //newItem.score = evaluateF1Score(&(cover), this->myclass);
                    newItem.score = (this->*evaluationFunction)(&(cover), this->myclass);            
                    newTmp.push(newItem);
                    ++iterate;

					int coverpos = 0;
					int coverneg = 0;
					getCoverage(&newItem, myclass, &coverpos, &coverneg);
					newItem.coverPos = coverpos;
					newItem.coverNeg = coverneg;
					//newItem.sig = sig;
                                        Rcpp::Rcout << "score: " << newItem.score << " pos: " << newItem.coverPos << " neg: " << newItem.coverNeg << " coverage: " << cover.count() << " mycandidate: ";
					for(int i = 0; i < newItem.rules.size(); ++i)
					{
						Node *mynode = newItem.rules[i]->noderef;
                                                Rcpp::Rcout << "id:#" << mynode->id << "# ";
					}

                                        Rcpp::Rcout << std::endl;
                    duplicitiesBITSET[getRuleBitset(&newItem)] = true;

                    //if(bestRule.score <= newItem.score)
                    //    bestRule = newItem;
                    double sig = evaluateSignificance(&newItem, myclass);
				//	if(sig > signTH)
				//	{
						atleastOneGenerated = true;
						if(newItem.score > bestRule.score)
							bestRule = newItem;
				//	}

                    //double entropy = evaluateEntropyScore(&(newItem), this->myclass);
                    //std::cout << entropy << std::endl;
                    //if(bestRuleEntropy.score > entropy)
                    //{
                    //    bestRuleEntropy = newItem;
                    //    bestRuleEntropy.score = entropy;
                    //}

                    //double acc = evaluateACCScore(&(newItem), this->myclass);
                    //std::cout << entropy << std::endl;
                    //if(bestRuleACC.score < acc)
                    //{
                    //    bestRuleACC = newItem;
                    //    bestRuleACC.score = acc;
                    //}
                }
            }
            filterComplexes(&newTmp, filterTH);
            ++iii;                        
            //std::cout << "i:" << iii << std::endl;
        }

        int beforeFilter = newTmp.size();
        filterComplexes(&newTmp, filterTH);
        newCandidates = newTmp;
        //std::cout << "total: " << beforeFilter << " runExhaustive: " << iterate  << " after filtering: "<< newCandidates.size() << " best score: " << bestRule.score << " rule: " << getPrintableFeature(&bestRule)  <<std::endl;
        //for(int ibest = 0; ibest < bestRule.rules.size(); ++ibest)
        //{
        //    std::cout << "id: " << bestRule.rules[ibest]->IDind << std::endl;
        //}
        Rcpp::Rcout << "RUN " << icycle+2 << " - size: " << iterate << " best score: " << bestRule.score << std::endl;
        Rcpp::Rcout << "Rule: " << getPrintableFeature(&bestRule) << std::endl << std::endl;
        clock_t end = std::clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        Rcpp::Rcout << "elapsed time: " << elapsed_secs << std::endl << std::endl;
        
       // if(!atleastOneGenerated)
	//	{
	//		return bestRule;
	//	}
        //std::cout << "Entropy " << " best score: " << bestRuleEntropy.score << " f1: " << evaluateF1Score(&bestRuleEntropy, this->myclass) << std::endl;
        //std::cout << "Rule: " << getPrintableFeature(&bestRuleEntropy) << std::endl << std::endl;

        //std::cout << "ACC " << " best score: " << bestRuleACC.score << " f1: " << evaluateF1Score(&bestRuleACC, this->myclass) << std::endl;
        //std::cout << "Rule: " << getPrintableFeature(&bestRuleACC) << std::endl << std::endl;


        //std::cout << std::endl << std::endl;
    }

/*    std::cout << "total score: ";
    while(!newCandidates.empty())
    {
        std::cout << newCandidates.top().score << " ";
        newCandidates.pop();
    }
*/
    return bestRule;
}

boost::dynamic_bitset<> BeamTopDown::getRuleBitset(newComplex *rule)
{
    boost::dynamic_bitset<> rulebitset(this->bottomFeatures->size(), 0); //rule->rules[0]->allChilds.size(),0);
    for(int ii = 0; ii < rule->rules.size(); ++ii)
    {
        rulebitset[rule->rules[ii]->IDind] = 1;
    }
    return rulebitset;
}

std::string BeamTopDown::createHashKey(std::vector<int> *identicalInt)
{
    std::sort(identicalInt->begin(), identicalInt->end());

    std::stringstream orderHash_stream;
    std::copy(identicalInt->begin(), identicalInt->end(), std::ostream_iterator<int>(orderHash_stream, " "));
    return orderHash_stream.str().c_str();
}

void BeamTopDown::filterComplexes(std::priority_queue<newComplex> *queue, int th)
{
    //remove elements
    while(queue->size() > th)
    {
        queue->pop();
    }
}

int BeamTopDown::runExhaustive(Rcpp::DataFrame *enrichPropositionTable)
{
	//https://teuder.gitbooks.io/introduction-to-rcpp/en/13_dataframe.html
    Rcpp::Rcout << "Starting exhaustive" << std::endl;
    Rcpp::NumericVector classColumn = (*enrichPropositionTable)["Class"];   //indexing start from 0
    int nsamples = classColumn.size();
    int ncol = (*enrichPropositionTable).length();
    boost::dynamic_bitset<> CLASS(nsamples,0);
    std::vector<boost::dynamic_bitset<> > groundFeatures(ncol - 1);

    //class numeric vector to bitset    
    Rcpp::Rcout << "class column size: " << nsamples << std::endl;
    for(int i = 0; i < nsamples; ++i)
    {
        if(classColumn[i] == 1)
            CLASS[i] = 1;
    }

    //feature numerical vectors to vector of bitsets
    for(int i = 0; i < ncol - 1; ++i)    //exclude class column
    {
        Rcpp::NumericVector featureVector = (*enrichPropositionTable)[i];
        boost::dynamic_bitset<> vector(nsamples, 0);
        for(int ibit = 0; ibit < featureVector.size(); ++ibit)
        {
            if(featureVector[ibit] == 1)
                vector[ibit] = 1;
        }
        groundFeatures[i] = vector;
    }


    std::vector<exhaustiveRule> F;
    //evaluate complex
    Rcpp::Rcout << "init score: ";
    for(int ieval = 0; ieval < groundFeatures.size(); ++ieval)
    {
	//std::cout << groundFeatures[ieval].count() << std::endl;
        double score = evaluateF1Score(&(groundFeatures[ieval]), &CLASS);	
        exhaustiveRule groud;
        groud.score = score;
        groud.coveredExample = groundFeatures[ieval];
        std::vector<boost::dynamic_bitset<> *> candidates;
        for(int iground = 0; iground < groundFeatures.size(); ++iground)
        {
            if(ieval != iground)
                candidates.push_back(&(groundFeatures[iground]));
        }

        groud.candidates = candidates;
        F.push_back(groud);
        Rcpp::Rcout << score << " ";
    }
    Rcpp::Rcout << std::endl;

    for(int icycle = 0; icycle < 4; ++icycle)
    {
        std::vector<exhaustiveRule> newF;
        Rcpp::Rcout << "cycle " << icycle << " score: ";
        for(int ground = 0; ground < F.size(); ++ground)
        {
            for(int add = 0; add < F[ground].candidates.size(); ++add)
            {
                boost::dynamic_bitset<> newRule = F[ground].coveredExample | (*F[ground].candidates[add]);
                double score = evaluateF1Score(&(newRule), &CLASS);
                exhaustiveRule newFeature;
                newFeature.score = score;
                newFeature.coveredExample = newRule;
                newFeature.candidates = F[ground].candidates;
                //erase current element
                newFeature.candidates.erase(newFeature.candidates.begin() + add);
                newF.push_back(newFeature);
                //std::cout << score << " ";
            }
        }
        Rcpp::Rcout << std::endl;
        Rcpp::Rcout << "newF size: " << newF.size() << std::endl;
        F = newF;
    }
}

newComplex BeamTopDown::runTopologyVersion(int filterTH, int ruleDepth, double signTH, int minLevel)
{
    //first step - initialize

    int iterate = 0;
    int totalIterates = 0;
    bool atleastOneGenerated = false;
    newComplex bestRule;
    bestRule.score = -1;


    clock_t begin = std::clock();
    std::priority_queue<newComplex> newCandidates;
    //bottom elements
    for(int ifeature = 0; ifeature < this->bottomFeatures->size(); ++ifeature)
    {
        newComplex newItem;
        //newItem.score = evaluateF1Score(&((*bottomFeatures)[ifeature].exampleCovered), this->myclass);
        //newItem.score = evaluateACCScore(&((*bottomFeatures)[ifeature].exampleCovered), this->myclass);
        //newItem.score = evaluateLaplaceScore(&((*bottomFeatures)[ifeature].exampleCovered), this->myclass);        
        newItem.score = (this->*evaluationFunction)(&((*bottomFeatures)[ifeature].exampleCovered), this->myclass);
        
        newItem.pathScore = newItem.score;
        newItem.rules.push_back(&((*bottomFeatures)[ifeature]));
        /*
        std::cout << "#####################" << std::endl; 
        std::cout << "ID: " << newItem.rules[0]->noderef->id  << std::endl;
        std::cout << "NAME: " << newItem.rules[0]->noderef->name  << std::endl;   
        std::cout << "score: " << newItem.score << std::endl;
	std::cout << "level: " << newItem.rules[0]->level << std::endl;
        std::cout << "#####################" << std::endl;
         */
        int coverpos = 0;
        int coverneg = 0;
        getCoverage(&newItem, this->myclass, &coverpos, &coverneg);
        newItem.coverPos = coverpos;
        newItem.coverNeg = coverneg;

        double sig = evaluateSignificance(&newItem, myclass);               
        newItem.sig = sig;
        
        newCandidates.push(newItem);
        //std::cout << "ssig: " << sig << std::endl;
        //if(sig > signTH)
        //{
			atleastOneGenerated = true;
            //if(bestRule.score <= newItem.score)
            //    bestRule = newItem;
        //}
        ++iterate;
    }

    if(this->verbose)
    {
        Rcpp::Rcout << "RUN 1" << " - size: " << iterate << " best score: " << bestRule.score << std::endl;
        Rcpp::Rcout << "Rule: " << getPrintableFeature(&bestRule) << std::endl << std::endl;
    }
    clock_t end = std::clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    if(this->verbose)
        Rcpp::Rcout << "elapsed time: " << elapsed_secs << std::endl << std::endl;


    filterComplexes(&newCandidates, filterTH);
    
	if(!atleastOneGenerated)
	{
		;//return bestRule;
	}

    boost::dynamic_bitset<> roots = Example::getRoots(this->bottomFeatures);
    if(this->verbose)
        Rcpp::Rcout << "roots: " << roots.count() << std::endl;
    std::vector<int> ontoRoots;
    size_t iroot = roots.find_first();
    while(iroot != boost::dynamic_bitset<>::npos)
    {
        ontoRoots.push_back(iroot);
        //next
        iroot = roots.find_next(iroot);
    }

    //length of rule
    for(int icycle = 0; icycle < ruleDepth; ++icycle)
    {
		atleastOneGenerated = false;
        std::priority_queue<newComplex> tmpNewCandidates;
        boost::unordered_map<boost::dynamic_bitset<>, bool> duplicitiesBITSET;        
        iterate = 0;
        int cf1 = 0;
        //for each element in queue
        while(!newCandidates.empty())
        {
            newComplex element2Expand = newCandidates.top();
            newCandidates.pop();


            //all childs and parents in rule, these nodes cannot be reached
            boost::dynamic_bitset<> allChilds = Example::getAllSpecifics(&element2Expand);
            boost::dynamic_bitset<> allParents = Example::getAllGenerals(&element2Expand);
            //std::cout << "allchilds: " << allChilds.count() << " allpar: " << allParents.count() << std::endl;

            boost::dynamic_bitset<> CLOSED_DAG = allChilds;//(allChilds.size(), 0);
            Example::addCurrNode2CLOSED(&element2Expand, &CLOSED_DAG);
            element2Expand.rules.resize(element2Expand.rules.size()+1);

            //for(int iii =0; iii < allChilds.size(); ++iii)
            //{
            //    if(allChilds[iii] == 1)
            //    {
            //        newComplex element2Expand2 = element2Expand;
            //        element2Expand2.rules[element2Expand2.rules.size()-1] = &((*bottomFeatures)[iii]);
            //        duplicitiesBITSET[getRuleBitset(&element2Expand2)] = true;
            //    }
            //}

            //size_t iallchild = allChilds.find_first();
            //while(iallchild != boost::dynamic_bitset<>::npos)
            //{
            //    newComplex element2Expand2 = element2Expand;
            //    element2Expand2.rules[element2Expand2.rules.size()-1] = &((*bottomFeatures)[iallchild]);
            //    duplicitiesBITSET[getRuleBitset(&element2Expand2)] = true;
            //    //next
            //    iallchild = allChilds.find_next(iallchild);
            //}

            //visit roots in a ontology
            //for(int ibit = 0; ibit < roots.size(); ++ibit)
            for(int iroot = 0; iroot < ontoRoots.size(); ++iroot)
            {
                //root
                int ibit = ontoRoots[iroot];
                if(CLOSED_DAG[ibit] == 0)
                {
                    //eval root
                    element2Expand.rules[element2Expand.rules.size()-1] = &((*bottomFeatures)[ibit]);
                    element2Expand.score = 0;//evaluateF1Score(&element2Expand, myclass);
                    element2Expand.pathScore = 0;//element2Expand.score;
                    //CLOSED[ibit] = 1;

                    //go down to DAG
                    std::list<newComplex> OPEN;
                    OPEN.push_back(element2Expand);
                    while(!OPEN.empty())
                    {
                        newComplex mycand = OPEN.front();
                        OPEN.pop_front();

                        //has not been visited
                        if(CLOSED_DAG[mycand.rules.back()->IDind] == 0)// && allChilds[mycand.rules.back()->IDind] == 0)
                        {
                            //is not parent
                            if(allParents[mycand.rules.back()->IDind] == 0 && duplicitiesBITSET.find(getRuleBitset(&mycand)) == duplicitiesBITSET.end())
                            {
                                double potencialScore;
                                double newscore;
                                //get f1 and potentail f1score
                                //evaluateF1ScorePotencialBoth(&mycand, myclass, &newscore, &potencialScore);
                                //evaluateACCScorePotencialBoth(&mycand, myclass, &newscore, &potencialScore);
                                
                                (this->*evaluationFunctionPotential)(&mycand, myclass, &newscore, &potencialScore);
                                
                                cf1++;
                                if(mycand.pathScore <= potencialScore && potencialScore >= bestRule.score)//1)//mycand.pathScore <= potencialScore) // TODO with LAPLACE!!!!! && potencialScore >= bestRule.score)
                                {
                                    //if(duplicitiesBITSET.find(getRuleBitset(&mycand)) == duplicitiesBITSET.end())
                                    //{

                                    //double scoreLaplace = evaluateLaplaceScore(&mycand, myclass);

                                    double sig = evaluateSignificance(&mycand, myclass);
                                    if(sig > signTH)
                                    {
                                        if(newscore > mycand.pathScore)
                                            mycand.pathScore = newscore;
                                    }
                                    mycand.score = newscore;//scoreLaplace;//newscore;
                                    
									//boost::dynamic_bitset<> cover(this->myclass->size());
									//cover.set();
									//bool hasColOntology = false;
									//for(int irule = 0; irule < mycand.rules.size(); ++irule)
									//{
										//std::cout << " b: " << mycand.rules[irule]->exampleCovered.count();
										//cover = cover & newItem.rules[irule]->exampleCovered;
									//	if(mycand.rules[irule]->noderef->ontoType != ROW_ONTOLOGY)
									//		hasColOntology = true;
									//}
                                    
                                    //if(!hasColOntology)
				//						mycand.score = 0;
                                    int coverpos = 0;
                                    int coverneg = 0;
                                    getCoverage(&mycand, myclass, &coverpos, &coverneg);
                                    mycand.coverPos = coverpos;
                                    mycand.coverNeg = coverneg;
                                    mycand.sig = sig;
									
                /*							Rcpp::Rcout << "score: " << mycand.score << " mycandidate: ";
									for(int i = 0; i < mycand.rules.size(); ++i)
									{
										Node *mynode = mycand.rules[i]->noderef;
										std::cout << "id:#" << mynode->id << "# ";
									}

									std::cout << std::endl;
          */                          tmpNewCandidates.push(mycand);
                                    duplicitiesBITSET[getRuleBitset(&mycand)] = true;
                                    iterate++;

                                    if(sig > signTH)
                                    {
                                        atleastOneGenerated = true;
                                        if(bestRule.score < mycand.score) //&& mycand.rules.back()->level >= minLevel)
					{
						bool isMinlevel = true;
						for(int ilevel = 0; ilevel < mycand.rules.size(); ++ilevel)
						{
							if(mycand.rules[ilevel]->level < minLevel)
							{
								isMinlevel = false;
								break;
							}
						}
						if(isMinlevel)
							bestRule = mycand;
					}
                                    }
                                    //}

                                }
                                else
                                {
                                    //has not potencial
                                    //CLOSED_DAG |= Example::completAllChildsLastRule(&mycand, bottomFeatures);
                                    //CLOSED_DAG[mycand.rules.back()->IDind] = 1;

                                    //addNonPotencialNodes2BITSET(&mycand, &duplicitiesBITSET);
                                    CLOSED_DAG |= Example::getAllSpecifics(&mycand);
                                    continue;
                                }
                            }

                            CLOSED_DAG[mycand.rules.back()->IDind] = 1;

                            //add siblings
                            size_t ichild = mycand.rules.back()->allChilds.find_first();
                            while(ichild != boost::dynamic_bitset<>::npos)
                            {
                                if(CLOSED_DAG[ichild] == 0)//allChilds[ichild] == 0 && CLOSED_DAG[ichild] == 0)
                                {
                                    newComplex newChild = mycand;
                                    newChild.rules[newChild.rules.size()-1] = &((*bottomFeatures)[ichild]);
                                    OPEN.push_back(newChild);
                                }
                                //next
                                ichild = mycand.rules.back()->allChilds.find_next(ichild);
                            }

                        }
                    }
                    filterComplexes(&tmpNewCandidates, filterTH);
                }
            }

        }
        filterComplexes(&tmpNewCandidates, filterTH);
        newCandidates = tmpNewCandidates;

        if(this->verbose)
        {
            Rcpp::Rcout << "RUN " << icycle+2 << " - size: " << iterate << " best score: " << bestRule.score << " count f1score: "  << cf1  << " cov stats: " << getPrintableCoverage(&bestRule, myclass) << std::endl;
            Rcpp::Rcout << "Rule: " << getPrintableFeature(&bestRule) << std::endl << std::endl;
        }
        clock_t end = std::clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

        if(this->verbose)
            Rcpp::Rcout << "elapsed time: " << elapsed_secs << std::endl << std::endl;
        
        if(!atleastOneGenerated)
        {
                ;//return bestRule;
        }
    }


    return bestRule;
}
std::vector<bottomFeature> BeamTopDown::makeComplex(newComplex *toExpand)
{
    ;
}

arma::mat BeamTopDown::removeCoveredExamples(newComplex actbest, arma::mat *new_armaData, boost::dynamic_bitset<> classBitMask)
{
    arma::mat actMatrix = *new_armaData;
    //std::cout << "BEFORE REMOVECOVERED: " << std::endl;
    //actMatrix.print();
    //arma::uvec indices = arma::find(actMatrix == 1);
    //indices.print();
    
    //std::cout << "MY POKUS!!!!!!" << std::endl;
    arma::uvec indicesPos = arma::find(actMatrix == 1);
    arma::umat tPos = arma::ind2sub(arma::size(actMatrix), indicesPos);
    
    arma::uvec indicesNeg = arma::find(actMatrix == 0);
    arma::umat tNeg = arma::ind2sub(arma::size(actMatrix), indicesNeg);
    
    int npositives = classBitMask.count();
    
    boost::dynamic_bitset<> rule = actbest.rules[0]->exampleCovered;
	for(int icomplex = 1; icomplex < actbest.rules.size(); ++icomplex)
	{
		rule &= actbest.rules[icomplex]->exampleCovered;
	}
	
	int myipos = 0;
	int myineg = 0;
	for(int ibit = 0; ibit < rule.size(); ++ibit)
	{
		if(rule[ibit] == 1)
		{			
			if(ibit < npositives)//indicesPos(indicesPos.size()-1))
			{
				//pos
				//actMatrix(ibit) = -1;
				int row = tPos.col(ibit)[0]; //row index
				int col = tPos.col(ibit)[1];  //col index
		//		std::cout << "positiv covered at: " << ibit << " position" << std::endl;
				actMatrix(row, col) = -1;
				++myipos;
			}
			else
			{
				//neg
				//actMatrix(ibit) = -1;
				int row = tNeg.col(ibit-npositives)[0]; //row index
				int col = tNeg.col(ibit-npositives)[1];  //col index
				actMatrix(row, col) = -1;
				++myineg;
			}
			//std::cout << ibit << " ";
		}
	}
    if(this->verbose)
    {
        Rcpp::Rcout << "pos: " << myipos << " neg: " << myineg << std::endl;
        Rcpp::Rcout << "cover: " << rule.count() << std::endl;
    }
    
    arma::uvec indicesPos1 = arma::find(actMatrix == 1);
    arma::uvec indicesNeg1 = arma::find(actMatrix == 0);
    arma::uvec indicesNeu = arma::find(actMatrix == -1);
    
    if(this->verbose)
        Rcpp::Rcout << "in new cycle pos: " << indicesPos1.size() << " neg: " << indicesNeg1.size() << " neut: " << indicesNeu.size() << std::endl << std::endl;
   // std::cout << "AFTER REMOVECOVERED: " << std::endl;
   // actMatrix.print();
    
    //actMatrix.print();
    /*
    std::cout << std::endl;
    for(int iElem = 0; iElem < indices2.size(); ++iElem)
        {
            std::cout << "i:" << iElem << " val:" << indices2(iElem) << std::endl;
        }
    
    
    
arma::umat tt       = arma::ind2sub(  arma::size(actMatrix), indices2 );
//t.print();
std::cout << std::endl;
for(int iElem = 0; iElem < tt.n_cols; ++iElem)
        {
            int irow = tt.col(iElem)[0]; //row index
            int icol = tt.col(iElem)[1];  //col index
            actMatrix(irow, icol) = -1; //covered
            std::cout << "row:" << irow << " col:" << icol << std::endl;
        }
        
        
        std::cout << std::endl;
        boost::dynamic_bitset<> rule = actbest.rules[0]->exampleCovered;
        for(int icomplex = 1; icomplex < actbest.rules.size(); ++icomplex)
        {
            rule &= actbest.rules[icomplex]->exampleCovered;
        }
        //std::cout << "cov:" << rule.count() << std::endl;
       // arma::uvec indices2(rule.size(), arma::fill::zeros);//(rule.size());


        for(int ibit = 0; ibit < rule.size(); ++ibit)
        {
            if(rule[ibit] == 1)
            {
				std::cout << ibit << " ";
                //indices[ibit] = ibit;
                //++tt;
            }
        }
        std::cout << std::endl;
*/
return actMatrix;
    /*
    if(actbest.rules.size() > 0)
    {
        //covered exampless
        boost::dynamic_bitset<> rule = actbest.rules[0]->exampleCovered;
        for(int icomplex = 1; icomplex < actbest.rules.size(); ++icomplex)
        {
            rule &= actbest.rules[icomplex]->exampleCovered;
        }
        std::cout << "cov:" << rule.count() << std::endl;
        arma::uvec indices(rule.size(), arma::fill::zeros);//(rule.size());


        int tt = 0;
        for(int ibit = 0; ibit < rule.size(); ++ibit)
        {
            if(rule[ibit] == 1)
            {
                indices[ibit] = ibit;
                ++tt;
            }
        }

        //indices.print();

        //std::cout << "tt: " << tt << std::endl;
        //size_t index = rule.find_first();
        //while(index != boost::dynamic_bitset::npos)
        //{
        //    indices[index] = 1;
        //    index = rule.find_next(index);
        //}
        //
        arma::umat indMatrix = arma::ind2sub(arma::size(actMatrix), indices);
        //indMatrix.print();
        std::cout << "size: " << arma::size(actMatrix) << std::endl;
        std::cout << "cov2: " << indMatrix.n_cols << std::endl;
        for(int iElem = 0; iElem < indMatrix.n_cols; ++iElem)
        {
            int irow = indMatrix.col(iElem)[0]; //row index
            int icol = indMatrix.col(iElem)[1];  //col index
            actMatrix(irow, icol) = -1; //covered
            //std::cout << "row:" << irow << " col:" << icol << std::endl;
        }
    }
    return actMatrix;
    * */
}

std::string BeamTopDown::printCoveredExamples(newComplex actbest, arma::mat *new_armaData, boost::dynamic_bitset<> classBitMask, Rcpp::CharacterVector rownames, Rcpp::CharacterVector colnames)
{
    arma::mat actMatrix = *new_armaData;
    arma::uvec indicesPos = arma::find(actMatrix == 1);
    arma::umat tPos = arma::ind2sub(arma::size(actMatrix), indicesPos);
    arma::uvec indicesNeg = arma::find(actMatrix == 0);
    arma::umat tNeg = arma::ind2sub(arma::size(actMatrix), indicesNeg);
    
    int npositives = classBitMask.count();
    
    boost::dynamic_bitset<> rule = actbest.rules[0]->exampleCovered;
	for(int icomplex = 1; icomplex < actbest.rules.size(); ++icomplex)
	{
		rule &= actbest.rules[icomplex]->exampleCovered;
	}
	
	int myipos = 0;
	int myineg = 0;
	std::string toprint = "  POSITIVE:\n";
	for(int ibit = 0; ibit < npositives; ++ibit)
	{
		if(rule[ibit] == 1)
		{			
			int row = tPos.col(ibit)[0]; //row index
			int col = tPos.col(ibit)[1];  //col index
			toprint += std::string("(") + Rcpp::as<std::string>(rownames(row)) + std::string(" - ");//std::string("(") + rownames(row) + std::string(" - ");
			toprint += Rcpp::as<std::string>(colnames(col)) +  std::string("), ");	
		}
	}
	
	toprint += "\n  NEGATIVE:\n";
	for(int ibit = npositives; ibit < rule.size(); ++ibit)
	{
		if(rule[ibit] == 1)
		{			
			int row = tNeg.col(ibit-npositives)[0]; //row index
			int col = tNeg.col(ibit-npositives)[1];  //col index
			toprint += std::string("(") + Rcpp::as<std::string>(rownames(row)) + std::string(" - ");
			toprint += Rcpp::as<std::string>(colnames(col)) + std::string("), ") ;	
		}
	}
	return toprint;
}

std::string BeamTopDown::printOnlyCoveredColExamples(newComplex actbest, arma::mat *new_armaData, boost::dynamic_bitset<> classBitMask, Rcpp::CharacterVector rownames, Rcpp::CharacterVector colnames)
{
    arma::mat actMatrix = *new_armaData;
    arma::uvec indicesPos = arma::find(actMatrix == 1);
    arma::umat tPos = arma::ind2sub(arma::size(actMatrix), indicesPos);
    arma::uvec indicesNeg = arma::find(actMatrix == 0);
    arma::umat tNeg = arma::ind2sub(arma::size(actMatrix), indicesNeg);
    
    int npositives = classBitMask.count();
    
    boost::dynamic_bitset<> rule = actbest.rules[0]->exampleCovered;
	for(int icomplex = 1; icomplex < actbest.rules.size(); ++icomplex)
	{
		rule &= actbest.rules[icomplex]->exampleCovered;
	}
	
	int myipos = 0;
	int myineg = 0;
	std::string toprint = "  POSITIVE:\n";
	for(int ibit = 0; ibit < npositives; ++ibit)
	{
		if(rule[ibit] == 1)
		{			
			//int row = tPos.col(ibit)[0]; //row index
			int col = tPos.col(ibit)[1];  //col index
			//toprint += std::string("(") + Rcpp::as<std::string>(rownames(row)) + std::string(" - ");//std::string("(") + rownames(row) + std::string(" - ");
			toprint += Rcpp::as<std::string>(colnames(col)) +  std::string(", ");
			ibit += rownames.size()-1;
		}
	}
	
	toprint += "\n  NEGATIVE:\n";
	for(int ibit = npositives; ibit < rule.size(); ++ibit)
	{
		if(rule[ibit] == 1)
		{			
			//int row = tNeg.col(ibit-npositives)[0]; //row index
			int col = tNeg.col(ibit-npositives)[1];  //col index
			//toprint += std::string("(") + Rcpp::as<std::string>(rownames(row)) + std::string(" - ");
			toprint += Rcpp::as<std::string>(colnames(col)) + std::string(", ") ;	
			ibit += rownames.size()-1;
		}
	}
	return toprint;
}


Rcpp::CharacterVector BeamTopDown::getOnlyCoveredRowExamples(newComplex actbest, arma::mat *new_armaData, Rcpp::CharacterVector rownames)
{
    Rcpp::CharacterVector coveredRules;
    arma::mat actMatrix = *new_armaData;
    arma::uvec indicesPos = arma::find(actMatrix == 1);
    arma::umat tPos = arma::ind2sub(arma::size(actMatrix), indicesPos);

    boost::dynamic_bitset<> rule = actbest.rules[0]->exampleCovered;
    for(int icomplex = 1; icomplex < actbest.rules.size(); ++icomplex)
    {
        rule &= actbest.rules[icomplex]->exampleCovered;
    }

    boost::unordered_map<int, bool> closed;
    for(int ibit = 0; ibit < indicesPos.n_rows; ++ibit)
    {
        int row = tPos.col(ibit)[0]; //row index
        int col = tPos.col(ibit)[1];

        if(rule[ibit] == 1)
        {
            if(closed.find(row) == closed.end())
            {
                coveredRules.push_back(Rcpp::as<std::string>(rownames(row)));
                closed[row] = true;
            }
        }
    }
    return coveredRules;
}

Rcpp::CharacterVector BeamTopDown::getOnlyCoveredColExamples(newComplex actbest, arma::mat *new_armaData, Rcpp::CharacterVector colnames)
{
    Rcpp::CharacterVector coveredRules;
    arma::mat actMatrix = *new_armaData;
    arma::uvec indicesPos = arma::find(actMatrix == 1);
    arma::umat tPos = arma::ind2sub(arma::size(actMatrix), indicesPos);

    boost::dynamic_bitset<> rule = actbest.rules[0]->exampleCovered;
    for(int icomplex = 1; icomplex < actbest.rules.size(); ++icomplex)
    {
        rule &= actbest.rules[icomplex]->exampleCovered;
    }

    boost::unordered_map<int, bool> closed;
    for(int ibit = 0; ibit < indicesPos.n_rows; ++ibit)
    {
        int row = tPos.col(ibit)[0]; //row index
        int col = tPos.col(ibit)[1];

        if(rule[ibit] == 1)
        {
            if(closed.find(col) == closed.end())
            {
                coveredRules.push_back(Rcpp::as<std::string>(colnames(col)));
                closed[col] = true;
            }
        }
    }
    return coveredRules;
}

Rcpp::CharacterVector BeamTopDown::getOnlyCoveredRowColExamples(newComplex actbest, arma::mat *new_armaData, Rcpp::CharacterVector rownames, Rcpp::CharacterVector colnames)
{
    Rcpp::CharacterVector coveredRules;
    arma::mat actMatrix = *new_armaData;
    arma::uvec indicesPos = arma::find(actMatrix == 1);
    arma::umat tPos = arma::ind2sub(arma::size(actMatrix), indicesPos);

    boost::dynamic_bitset<> rule = actbest.rules[0]->exampleCovered;
    for(int icomplex = 1; icomplex < actbest.rules.size(); ++icomplex)
    {
        rule &= actbest.rules[icomplex]->exampleCovered;
    }

    for(int ibit = 0; ibit < indicesPos.n_rows; ++ibit)
    {
        int row = tPos.col(ibit)[0]; //row index
        int col = tPos.col(ibit)[1];

        if(rule[ibit] == 1)
        {
                std::string desc = Rcpp::as<std::string>(rownames(row));
                desc += ",";
                desc += colnames(col);
                coveredRules.push_back(desc);
        }
    }
    return coveredRules;
}

void BeamTopDown::addNonPotencialNodes2BITSET(newComplex *mycand, boost::unordered_map<boost::dynamic_bitset<>, bool> *duplicitiesBITSET)
{
    boost::dynamic_bitset<> allSpecific = Example::getAllSpecifics(mycand);
    newComplex firstNode = *mycand;
    (*duplicitiesBITSET)[getRuleBitset(&firstNode)] = true;

    size_t ichild = allSpecific.find_first();
    while(ichild != boost::dynamic_bitset<>::npos)
    {
        firstNode.rules[firstNode.rules.size()-1] = &((*bottomFeatures)[ichild]);
        (*duplicitiesBITSET)[getRuleBitset(&firstNode)] = true;
        //next
        ichild = allSpecific.find_next(ichild);
    }

    //for(int i = 0; i < allSpecific.size(); ++i)
    //{
    //    if(allSpecific[i] == 1)
    //    {
    //        firstNode.rules[firstNode.rules.size()-1] = &((*bottomFeatures)[i]);
    //        (*duplicitiesBITSET)[getRuleBitset(&firstNode)] = true;
    //    }
    //}
}


std::string BeamTopDown::getPrintableFeature(newComplex *toPrint)
{
    std::string rule = "";
    if(toPrint->rules.size() > 0)
    {
        if(toPrint->rules[0]->posFeature)
            rule += toPrint->rules[0]->noderef->name + std::string("(") + boost::lexical_cast<std::string>(toPrint->rules[0]->noderef->idNum) + std::string(")");
        else
            rule += "!" + toPrint->rules[0]->noderef->name + std::string("(") + boost::lexical_cast<std::string>(toPrint->rules[0]->noderef->idNum) + std::string(")");
        for(int irule = 1; irule < toPrint->rules.size(); ++irule)
        {
            //negative
            if(toPrint->rules[irule]->posFeature)
                rule += "\n" + toPrint->rules[irule]->noderef->name + std::string("(") + boost::lexical_cast<std::string>(toPrint->rules[irule]->noderef->idNum) + std::string(")");
            else
                rule += "\n!" + toPrint->rules[irule]->noderef->name + std::string("(") + boost::lexical_cast<std::string>(toPrint->rules[irule]->noderef->idNum) + std::string(")");
        }
    }
    return rule;
}

std::string BeamTopDown::getPrintableCoverage(newComplex *toPrint, boost::dynamic_bitset<> *classVector)
{
    std::string rulestring = "";
    if(toPrint->rules.size() > 0)
    {
        boost::dynamic_bitset<> rule = toPrint->rules[0]->exampleCovered;
        for(int icomplex = 0; icomplex < toPrint->rules.size(); ++icomplex)
        {
            rule &= toPrint->rules[icomplex]->exampleCovered;
        }
        int positive = (rule & *classVector).count();
        int negative = (rule & ~(*classVector)).count();
        rulestring += std::string("[") + boost::lexical_cast<std::string>(positive) + std::string(", ") + boost::lexical_cast<std::string>(negative) + std::string("]");
    }
    return rulestring;
}

void BeamTopDown::getCoverage(newComplex *toPrint, boost::dynamic_bitset<> *classVector, int *npos, int *nneg)
{
    if(toPrint->rules.size() > 0)
    {
        boost::dynamic_bitset<> rule = toPrint->rules[0]->exampleCovered;
        for(int icomplex = 1; icomplex < toPrint->rules.size(); ++icomplex)
        {
            rule &= toPrint->rules[icomplex]->exampleCovered;
        }
        *npos = (rule & *classVector).count();
        *nneg = (rule & ~(*classVector)).count();       
    }
}

double BeamTopDown::evaluateF1Score(boost::dynamic_bitset<> *vector, boost::dynamic_bitset<> *classVector)
{
    int tp = computeTP(vector, classVector);
    if(tp == 0)
	return 0;
    else
    {
	double denominator = 2*tp + computeFP(vector, classVector) + computeFN(vector, classVector);
	return (2*tp)/denominator;
    }
    return (2*(computeTP(vector, classVector)))/(double)(2*(computeTP(vector, classVector)) + computeFP(vector, classVector) + (computeFN(vector, classVector)));
}

double BeamTopDown::evaluateF1ScorePotencial(boost::dynamic_bitset<> *vector, boost::dynamic_bitset<> *classVector)
{
    int tp = computeTP(vector, classVector);
    if(tp == 0)
        return 0;
    else
    {
        double denominator = 2*tp + computeFN(vector, classVector);
        return (2*tp)/denominator;
    }
    return (2*(computeTP(vector, classVector)))/(double)(2*(computeTP(vector, classVector)) + computeFP(vector, classVector) + (computeFN(vector, classVector)));
}

void BeamTopDown::evaluateF1ScorePotencialBoth(boost::dynamic_bitset<> *vector, boost::dynamic_bitset<> *classVector, double *f1score, double *potential)
{
    int tp = computeTP(vector, classVector);
    if(tp == 0)
    {
        *f1score = 0;
        *potential = 0;
    }
    else
    {
        double denominatorPot = 2*tp + computeFN(vector, classVector);
        double denominatorF1 = denominatorPot + computeFP(vector, classVector);
        *f1score = (2*tp)/denominatorF1;
        *potential = (2*tp)/denominatorPot;
    }
}

void BeamTopDown::evaluateACCScorePotencialBoth(boost::dynamic_bitset<> *vector, boost::dynamic_bitset<> *classVector, double *accscore, double *potential)
{
    int tp = computeTP(vector, classVector);
    int tn = computeTN(vector, classVector);
    if(tp + tn == 0)
    {
        *accscore = 0;
        *potential = 0;
    }
    else
    {
        int fn = computeFN(vector, classVector);
        int fp = computeFP(vector, classVector);
        double denom = tp + tn + fn + fp;
        //double denominatorPot = tp + tn + computeFN(vector, classVector);
        //double denominatorACC = denominatorPot + computeFP(vector, classVector);
        int numerator = tp + tn;
        int numeratorPot = numerator + fp;
        //*accscore = (tp+tn)/denominatorACC;
        //*potential = (tp+tn)/denominatorPot;

        *accscore = numerator/denom;
        *potential = numeratorPot/denom;
    }
}

double BeamTopDown::evaluateF1Score(newComplex *complex, boost::dynamic_bitset<> *classVector)
{
    if(complex->rules.size() > 0)
    {
        boost::dynamic_bitset<> rule = complex->rules[0]->exampleCovered;
        for(int icomplex = 0; icomplex < complex->rules.size(); ++icomplex)
        {
            rule &= complex->rules[icomplex]->exampleCovered;
        }
        return evaluateF1Score(&rule, classVector);
    }
    else
        return 0;
}

double BeamTopDown::evaluateF1ScorePotencial(newComplex *complex, boost::dynamic_bitset<> *classVector)
{
    if(complex->rules.size() > 0)
    {
        boost::dynamic_bitset<> rule = complex->rules[0]->exampleCovered;
        for(int icomplex = 0; icomplex < complex->rules.size(); ++icomplex)
        {
            rule &= complex->rules[icomplex]->exampleCovered;
        }
        return evaluateF1ScorePotencial(&rule, classVector);
    }
    else
        return 0;
}

void BeamTopDown::evaluateF1ScorePotencialBoth(newComplex *complex, boost::dynamic_bitset<> *classVector, double *f1score, double *potencial)
{
    if(complex->rules.size() > 0)
    {
        boost::dynamic_bitset<> rule = complex->rules[0]->exampleCovered;
        for(int icomplex = 1; icomplex < complex->rules.size(); ++icomplex)
        {
            rule &= complex->rules[icomplex]->exampleCovered;
        }
        //actualize f1 and potentail f1score
        evaluateF1ScorePotencialBoth(&rule, classVector, f1score, potencial);
    }
    else
    {
        *f1score = 0;
        *potencial = 0;
    }
}

void BeamTopDown::evaluateACCScorePotencialBoth(newComplex *complex, boost::dynamic_bitset<> *classVector, double *accscore, double *potencial)
{
    if(complex->rules.size() > 0)
    {
        boost::dynamic_bitset<> rule = complex->rules[0]->exampleCovered;
        for(int icomplex = 1; icomplex < complex->rules.size(); ++icomplex)
        {
            rule &= complex->rules[icomplex]->exampleCovered;
        }
        //actualize f1 and potentail f1score
        evaluateACCScorePotencialBoth(&rule, classVector, accscore, potencial);
    }
    else
    {
        *accscore = 0;
        *potencial = 0;
    }
}

void BeamTopDown::evaluateAUCScorePotencialBoth(boost::dynamic_bitset<> *vector, boost::dynamic_bitset<> *classVector, double *score, double *potential)
{
    int tp = computeTP(vector, classVector);
    int tn = computeTN(vector, classVector);
    int fp = computeFP(vector, classVector);
    int fn = computeFN(vector, classVector);

    double tpr = 0;
    double fpr = 0;
    if(tp > 0 || fn > 0)
		tpr = tp/(double)(tp+fn);
	if(fp > 0 || tn > 0)
		fpr = fp/(double)(fp+tn);

    *score = (fpr*tpr)/2 + (1-fpr)*tpr + ((1-fpr)*(1-tpr))/2;
	*potential = tpr + (1-tpr)/2;	//fpr is zero    
}


void BeamTopDown::evaluateAUCScorePotencialBoth(newComplex *complex, boost::dynamic_bitset<> *classVector, double *score, double *potencial)
{
    if(complex->rules.size() > 0)
    {
        boost::dynamic_bitset<> rule = complex->rules[0]->exampleCovered;
        for(int icomplex = 1; icomplex < complex->rules.size(); ++icomplex)
        {
            rule &= complex->rules[icomplex]->exampleCovered;
        }
        //actualize f1 and potentail f1score
        evaluateAUCScorePotencialBoth(&rule, classVector, score, potencial);
    }
    else
    {
        *score = 0;
        *potencial = 0;
    }
}


double BeamTopDown::evaluateEntropyScore(boost::dynamic_bitset<> *vector, boost::dynamic_bitset<> *classVector)
{
    int tp = computeTP(vector, classVector);
    int tn = computeTN(vector, classVector);
    double laplace = 0.00000001;

    double pos = tp / (double)classVector->count();
    pos += laplace;
    double neg = tn / (double)(~(*classVector)).count();
    neg += laplace;
    return -((pos)*log2(pos) + (neg)*log2(neg));
}

double BeamTopDown::evaluateEntropyScore(newComplex *complex, boost::dynamic_bitset<> *classVector)
{
    if(complex->rules.size() > 0)
    {
        boost::dynamic_bitset<> rule = complex->rules[0]->exampleCovered;
        for(int icomplex = 0; icomplex < complex->rules.size(); ++icomplex)
        {
            rule &= complex->rules[icomplex]->exampleCovered;
        }
        return evaluateEntropyScore(&rule, classVector);
    }
    else
        return 0;
}

double BeamTopDown::evaluateACCScore(boost::dynamic_bitset<> *vector, boost::dynamic_bitset<> *classVector)
{
    int tp = computeTP(vector, classVector);
    int tn = computeTN(vector, classVector);
    int fp = computeFP(vector, classVector);
    int fn = computeFN(vector, classVector);

    return (tp+tn)/(double)(tp+tn+fp+fn);
}

double BeamTopDown::evaluateAUCcore(boost::dynamic_bitset<> *vector, boost::dynamic_bitset<> *classVector)
{
    int tp = computeTP(vector, classVector);
    int tn = computeTN(vector, classVector);
    int fp = computeFP(vector, classVector);
    int fn = computeFN(vector, classVector);

	double tpr = 0;
	double fpr = 0;
	if(tp > 0 || fn > 0)
		tpr = tp/(double)(tp+fn);
	if(fp > 0 || tn > 0)
		fpr = fp/(double)(fp+tn);

    return (fpr*tpr)/2 + (1-fpr)*tpr + ((1-fpr)*(1-tpr))/2;
}

double BeamTopDown::evaluateLaplaceScore(boost::dynamic_bitset<> *vector, boost::dynamic_bitset<> *classVector)
{
    int tp = computeTP(vector, classVector);
    int tpfp = vector->count();
    return (tp+1)/(double)(tpfp+2);
}

double BeamTopDown::evaluateLaplaceScore(newComplex *complex, boost::dynamic_bitset<> *classVector)
{
    if(complex->rules.size() > 0)
    {
        boost::dynamic_bitset<> rule = complex->rules[0]->exampleCovered;
        for(int icomplex = 0; icomplex < complex->rules.size(); ++icomplex)
        {
            rule &= complex->rules[icomplex]->exampleCovered;
        }
        return evaluateLaplaceScore(&rule, classVector);
    }
    else
        return 0;
}

double BeamTopDown::evaluateACCScore(newComplex *complex, boost::dynamic_bitset<> *classVector)
{
    if(complex->rules.size() > 0)
    {
        boost::dynamic_bitset<> rule = complex->rules[0]->exampleCovered;
        for(int icomplex = 0; icomplex < complex->rules.size(); ++icomplex)
        {
            rule &= complex->rules[icomplex]->exampleCovered;
        }
        return evaluateACCScore(&rule, classVector);
    }
    else
        return 0;
}

double BeamTopDown::evaluateSignificance(newComplex *complex, boost::dynamic_bitset<> *classVector)
{
    if(complex->rules.size() > 0)
    {
        boost::dynamic_bitset<> rule = complex->rules[0]->exampleCovered;
        for(int icomplex = 1; icomplex < complex->rules.size(); ++icomplex)
        {
            rule &= complex->rules[icomplex]->exampleCovered;
        }
        double tp = computeTP(&rule, classVector);
        //double tn = computeTN(&rule, classVector);
        double fp = computeFP(&rule, classVector);
        double posExamples = classVector->count();
        double negExamples = (~(*classVector)).count();
        //double posCovered = rule.count();
        //double negCovered = (~rule).count();
        //double n = classVector->size();

        double res1 = tp * log2((tp/(tp+fp))/(posExamples/(posExamples+negExamples)));
        double res2 = fp * log2((fp/(tp+fp))/(negExamples/(posExamples+negExamples)));
        return 2*(res1+res2);

        /*if(n != 0)
        {

            double resPos = 0;
            double resNeg = 0;
            if(posCovered != 0 && posExamples != 0 && tp != 0)
                resPos = tp*log2((double)posCovered*((double)posCovered/n));//tp/(posExamples*(posCovered/(n)));//resPos = tp*log2(tp/(posExamples*(posCovered/(n))));
            if(negCovered != 0 && negExamples != 0 && tn != 0)
                resNeg = tn*log2((double)negCovered*((double)negCovered/n));//tn/(negExamples*(negCovered/(n)));//resNeg = tn*log2(tn/(negExamples*(negCovered/(n))));
            if(resPos == 0 && resNeg == 0)
                return 1;
            else
                return 2*(resPos  + resNeg);
        }
        */
    }
    return 1;
}

int BeamTopDown::computeTP(boost::dynamic_bitset<> *vector, boost::dynamic_bitset<> *classVector)
{
    return (*vector & *classVector).count();
}

int BeamTopDown::computeFP(boost::dynamic_bitset<> *vector, boost::dynamic_bitset<> *classVector)
{
    return (*vector & ~(*classVector)).count();
}

int BeamTopDown::computeFN(boost::dynamic_bitset<> *vector, boost::dynamic_bitset<> *classVector)
{
    return (~(*vector) & *classVector).count();
}

int BeamTopDown::computeTN(boost::dynamic_bitset<> *vector, boost::dynamic_bitset<> *classVector)
{
    return (~(*vector) & ~(*classVector)).count();
}

int BeamTopDown::run()
{

    std::priority_queue<conjunction_max> max_heap;
    std::priority_queue<conjunction_max> max_heap_TOTAL;
    std::priority_queue<conjunction_min> min_heap;
    std::priority_queue<conjunction_max> OPEN_heap;

    conjunction_max ABSL_BEST;
    ABSL_BEST.object = 0;

    int BEAM_LIMIT = 300;

    //CLOSED for each ontology
    std::vector<boost::dynamic_bitset<> > CLOSED;
    for(int i = 0; i < refOntologies->size(); ++i)
    {
        boost::dynamic_bitset<> ontoBit((*refOntologies)[i]->getOntologyParser()->getTermsNumber(), 0);
        CLOSED.push_back(ontoBit);
    }


    //INIT PHASE

    for(int iOnto = 0; iOnto < this->refOntologies->size(); ++iOnto)
    {
        std::vector<Node*> *roots = (*refOntologies)[iOnto]->getOntologyParser()->getRoots();
        //from bottom up to TOP DOWN
        for(int i = 0; i < roots->size(); ++i)
        {

            (*roots)[i]= (*roots)[i]->onto_ref->getOntologyParser()->getNodeTopDown(&((*roots)[i]->id));
        }

        for(int iroot = 0; iroot < roots->size(); ++iroot)
        {
            //initialization
            conjunction_max newItem;
            newItem.node2Specified = (*roots)[iroot];   //set current node!
            newItem.bitset = std::vector<boost::dynamic_bitset<> >(this->refOntologies->size());
            newItem.ontoID = std::vector<int>(this->refOntologies->size(),0);
            std::vector<boost::dynamic_bitset<> > actDiscovered(this->refOntologies->size());
            for(int ihypothesis = 0; ihypothesis < this->refOntologies->size(); ++ihypothesis)
            {
                newItem.ontoID[ihypothesis] = ihypothesis;
                boost::dynamic_bitset<> bitset((*refOntologies)[ihypothesis]->getOntologyParser()->getTermsNumber(),0);
                newItem.bitset[ihypothesis] = bitset;

                actDiscovered[ihypothesis] = ~(*this->enrichItems)[ihypothesis];
            }
            newItem.discovered = actDiscovered;
            newItem.bitset[iOnto][(*roots)[iroot]->idNum]  = 1; //set root in the current ontology
            newItem.intNode2Specified = iOnto;
            newItem.object = this->evaluateF1Score(&newItem);//newItem.bitset[ontoID], &ontoID);
            newItem.bestPathScore = newItem.object;
            Rcpp::Rcout << "root obj: " << newItem.object  << " ONTO: " << (*roots)[iroot]->onto_ref->getName() << " P: " << newItem.P << " N: " << newItem.N << std::endl;
            CLOSED[iOnto][(*roots)[iroot]->idNum] = 1;//set closed



            //max_heap.push(newItem);
            OPEN_heap.push(newItem);
            for(int i = 0; i < ITERATE_LIMIT; ++i)
            {
                if(OPEN_heap.empty())
                        break;

                generateNewConjunctions(&OPEN_heap, &max_heap, &CLOSED);

            }
            if(!QUEUE_UNLIMITED)
                checkQueueMaxLimit(&max_heap, &min_heap, 11);

            if(ABSL_BEST < max_heap.top())
                ABSL_BEST = max_heap.top();

            while(!max_heap.empty())
            {
                conjunction_max best = max_heap.top();
                max_heap.pop();
                max_heap_TOTAL.push(best);
                //topBest.push_back(best);
                Rcpp::Rcout << "onto: " << best.node2Specified->onto_ref->getName() << " score: " << best.object << " P:" << best.P << " N:" << best.N  << " hypothesis: ";
                best.node2Specified->onto_ref->printSemanticPattern(&(best.bitset[best.intNode2Specified]));

                for(int io = 0; io < refOntologies->size(); ++io)
                {
                    if(!best.bitset[io].none())
                    {
                        Rcpp::Rcout << (*refOntologies)[io]->getName() << "_";
                        (*refOntologies)[io]->printSemanticPattern(&(best.bitset[io]));
                    }
                }

                Rcpp::Rcout << std::endl << "--------------------------" << std::endl;
            }

            while(!OPEN_heap.empty())
            {
                OPEN_heap.pop();
            }
        }
    }

    Rcpp::Rcout << "STOP INITIALIZING"  << std::endl;

    for(int iicycle = 0; iicycle < 300; ++iicycle)
    {
        conjunction_max actbest = max_heap_TOTAL.top();
        max_heap_TOTAL.pop();
        //GET BEST and EVALUATE THEM
        for(int iOnto = 0; iOnto < this->refOntologies->size(); ++iOnto)
        {
            std::vector<Node*> *roots = (*refOntologies)[iOnto]->getOntologyParser()->getRoots();
            //from bottom up to TOP DOWN
            for(int i = 0; i < roots->size(); ++i)
            {

                (*roots)[i]= (*roots)[i]->onto_ref->getOntologyParser()->getNodeTopDown(&((*roots)[i]->id));
            }

            for(int iroot = 0; iroot < roots->size(); ++iroot)
            {
                //initialization
                //conjunction_max newItem;

                for(int i = 0; i < refOntologies->size(); ++i)
                {
                    boost::dynamic_bitset<> ontoBit((*refOntologies)[i]->getOntologyParser()->getTermsNumber(), 0);
                    CLOSED[i] = ontoBit;
                }

                conjunction_max newItem = actbest;

                newItem.node2Specified = (*roots)[iroot];   //set current node!
                newItem.bitset[iOnto][(*roots)[iroot]->idNum]  = 1; //set root in the current ontology
                newItem.intNode2Specified = iOnto;
                newItem.object = this->evaluateF1Score(&newItem);//newItem.bitset[ontoID], &ontoID);
                newItem.bestPathScore = newItem.object;
                Rcpp::Rcout << "new - root obj (" << (*roots)[iroot]->onto_ref->getName() << ") : " << newItem.object << " P: " << newItem.P << " N: " << newItem.N << std::endl;


                Rcpp::Rcout << "to discover- onto: " << newItem.node2Specified->onto_ref->getName()  <<  " score: " << newItem.object << " P:" << newItem.P << " N:" << newItem.N  << " hypothesis: ";
                newItem.node2Specified->onto_ref->printSemanticPattern(&(newItem.bitset[newItem.intNode2Specified]));
                Rcpp::Rcout << std::endl;


                CLOSED[iOnto][(*roots)[iroot]->idNum] = 1;//set closed

                OPEN_heap.push(newItem);
                for(int i = 0; i < ITERATE_LIMIT; ++i)
                {
                    if(OPEN_heap.empty())
                            break;

                    generateNewConjunctions(&OPEN_heap, &max_heap, &CLOSED);

                }
                if(!QUEUE_UNLIMITED)
                    checkQueueMaxLimit(&max_heap, &min_heap, 11);

                if(ABSL_BEST < max_heap.top())
                    ABSL_BEST = max_heap.top();

                while(!max_heap.empty())
                {
                    conjunction_max best = max_heap.top();
                    max_heap.pop();
                    max_heap_TOTAL.push(best);
                    //topBest.push_back(best);
                    Rcpp::Rcout << "onto: " << best.node2Specified->onto_ref->getName() << " score: " << best.object << " P:" << best.P << " N:" << best.N  << " hypothesis: ";
                    best.node2Specified->onto_ref->printSemanticPattern(&(best.bitset[best.intNode2Specified]));

                    for(int io = 0; io < refOntologies->size(); ++io)
                    {
                        if(!best.bitset[io].none())
                        {
                            Rcpp::Rcout << (*refOntologies)[io]->getName() << "_";
                            (*refOntologies)[io]->printSemanticPattern(&(best.bitset[io]));
                        }
                    }

                    Rcpp::Rcout << std::endl << "--------------------------" << std::endl;
                }

                while(!OPEN_heap.empty())
                {
                    OPEN_heap.pop();
                }



            }
        }
    }

    Rcpp::Rcout << " --- TOTAL ---" << std::endl;
    while(!max_heap_TOTAL.empty())
    {
        conjunction_max best = max_heap_TOTAL.top();
        max_heap_TOTAL.pop();
        Rcpp::Rcout << "onto: " << best.node2Specified->onto_ref->getName() << " score: " << best.object << " P:" << best.P << " N:" << best.N  << " hypothesis: ";
        best.node2Specified->onto_ref->printSemanticPattern(&(best.bitset[best.intNode2Specified]));

        for(int io = 0; io < refOntologies->size(); ++io)
        {
            if(!best.bitset[io].none())
            {
                Rcpp::Rcout << (*refOntologies)[io]->getName() << "_";
                (*refOntologies)[io]->printSemanticPattern(&(best.bitset[io]));
            }
        }

        Rcpp::Rcout << std::endl; //<< "--------------------------" << std::endl;
    }

    Rcpp::Rcout << std::endl << "BEST OF BEST: ";
    Rcpp::Rcout << "onto: " << ABSL_BEST.node2Specified->onto_ref->getName() << " score: " << ABSL_BEST.object << " P:" << ABSL_BEST.P << " N:" << ABSL_BEST.N  << " hypothesis: ";
    //ABSL_BEST.node2Specified->onto_ref->printSemanticPattern(&(ABSL_BEST.bitset[ABSL_BEST.intNode2Specified]));

    for(int io = 0; io < refOntologies->size(); ++io)
    {
        if(!ABSL_BEST.bitset[io].none())
        {
            Rcpp::Rcout << (*refOntologies)[io]->getName() << "_";
            (*refOntologies)[io]->printSemanticPattern(&(ABSL_BEST.bitset[io]));
        }
    }


    return 1;
    //END INIT PHASE






    conjunction_max bestConjunction;
    //first step
    initPriorityQueue(&OPEN_heap, &max_heap, &CLOSED);
    conjunction_max tmpBest = OPEN_heap.top();
    OPEN_heap.pop();

    max_heap_TOTAL = OPEN_heap;

    while(!OPEN_heap.empty())
    {
        OPEN_heap.pop();
    }
    OPEN_heap.push(tmpBest);

    for(int icycles = 0; icycles < 2; ++icycles)
    {
        conjunction_max bestConjunctionTMP = OPEN_heap.top();
        Rcpp::Rcout << "to discover- onto: " << bestConjunctionTMP.node2Specified->onto_ref->getName()  <<  " score: " << bestConjunctionTMP.object << " P:" << bestConjunctionTMP.P << " N:" << bestConjunctionTMP.N  << " hypothesis: ";
        bestConjunctionTMP.node2Specified->onto_ref->printSemanticPattern(&(bestConjunctionTMP.bitset[bestConjunctionTMP.intNode2Specified]));
        Rcpp::Rcout << std::endl;
        for(int i = 0; i < ITERATE_LIMIT; ++i)
        {
            if(OPEN_heap.empty())
                    break;

            generateNewConjunctions(&OPEN_heap, &max_heap, &CLOSED);
            //if(!QUEUE_UNLIMITED)
            //        checkQueueMaxLimit(&max_heap, &min_heap);
            bestConjunction = max_heap.top();

        }
        Rcpp::Rcout << "Best generated item: " << std::endl;
        int pos = 1;
        std::vector<conjunction_max> topBest;
        //check limit
        if(!QUEUE_UNLIMITED)
            checkQueueMaxLimit(&max_heap, &min_heap, 2);

        while(!max_heap.empty())
        {
            conjunction_max best = max_heap.top();
            max_heap.pop();
            max_heap_TOTAL.push(best);
            //topBest.push_back(best);
            Rcpp::Rcout << "onto: " << best.node2Specified->onto_ref->getName() << " score: " << best.object << " P:" << best.P << " N:" << best.N  << " hypothesis: ";
            best.node2Specified->onto_ref->printSemanticPattern(&(best.bitset[best.intNode2Specified]));

            for(int io = 0; io < refOntologies->size(); ++io)
            {
                if(!best.bitset[io].none())
                {
                    Rcpp::Rcout << (*refOntologies)[io]->getName() << "_";
                    (*refOntologies)[io]->printSemanticPattern(&(best.bitset[io]));
                }
            }

            Rcpp::Rcout << std::endl;
            pos++;
        }
        Rcpp::Rcout << std::endl;
      /*  for(int ib = 0; ib < topBest.size(); ++ib)
        {
            max_heap.push(topBest[ib]);
        }
        */
        while(!OPEN_heap.empty())
        {
            OPEN_heap.pop();
        }

        for(int iclosed = 0; iclosed < CLOSED.size(); ++iclosed)
        {
            CLOSED[iclosed].reset();
        }


        conjunction_max actbest = max_heap_TOTAL.top();
        max_heap_TOTAL.pop();

        Rcpp::Rcout << "Best of ALL: ";
        Rcpp::Rcout << "onto: " << actbest.node2Specified->onto_ref->getName() << " score: " << actbest.object << " P:" << actbest.P << " N:" << actbest.N  << " hypothesis: ";
        actbest.node2Specified->onto_ref->printSemanticPattern(&(actbest.bitset[actbest.intNode2Specified]));
        Rcpp::Rcout << std::endl;

        for(int iOnto = 0; iOnto < this->refOntologies->size(); ++iOnto)
        {
            std::vector<Node*> *roots = (*refOntologies)[iOnto]->getOntologyParser()->getRoots();
            //from bottom up to TOP DOWN
            for(int i = 0; i < roots->size(); ++i)
            {

                (*roots)[i]= (*roots)[i]->onto_ref->getOntologyParser()->getNodeTopDown(&((*roots)[i]->id));
            }

            for(int iroot = 0; iroot < roots->size(); ++iroot)
            {
                //initialization
                //conjunction_max newItem;
                conjunction_max newItem = actbest;
                newItem.node2Specified = (*roots)[iroot];   //set current node!
                //newItem.bitset = std::vector<boost::dynamic_bitset<> >(this->refOntologies->size());
                //newItem.ontoID = std::vector<int>(this->refOntologies->size(),0);
                //for(int ihypothesis = 0; ihypothesis < this->refOntologies->size(); ++ihypothesis)
                //{
                //    newItem.ontoID[ihypothesis] = ihypothesis;
                //    boost::dynamic_bitset<> bitset((*refOntologies)[ihypothesis]->getOntologyParser()->getTermsNumber(),0);
                //    newItem.bitset[ihypothesis] = bitset;
                //}
                newItem.bitset[iOnto][(*roots)[iroot]->idNum]  = 1; //set root in the current ontology
                newItem.intNode2Specified = iOnto;
                newItem.object = this->evaluateF1Score(&newItem);//newItem.bitset[ontoID], &ontoID);
                newItem.bestPathScore = newItem.object;
                Rcpp::Rcout << "new - root obj (" << (*roots)[iroot]->onto_ref->getName() << ") : " << newItem.object << " P: " << newItem.P << " N: " << newItem.N << std::endl;
                //boost::dynamic_bitset<> *bitset = new boost::dynamic_bitset<>((*refOntologies)[iOnto]->getOntologyParser()->getTermsNumber(),0);
                //(*bitset)[(*roots)[iroot]->idNum] = 1;    //set
                CLOSED[iOnto][(*roots)[iroot]->idNum] = 1;//set closed

                //generate but cannot be in a hypothesis
                //if(newItem.discovered[newItem.intNode2Specified][(*roots)[iroot]->idNum] != 1)
                //{
                //    max_heap.push(newItem);
                    max_heap_TOTAL.push(newItem);
                //}


                //OPEN_heap.push(newItem);
            }
        }

        //check limit
        if(!QUEUE_UNLIMITED)
            checkQueueMaxLimit(&max_heap_TOTAL, &min_heap, 100);


        Rcpp::Rcout << "CURRENT TOTAL SCORE:" << std::endl;
        topBest.clear();

        while(!max_heap_TOTAL.empty())
        {
            conjunction_max best = max_heap_TOTAL.top();
            max_heap_TOTAL.pop();
            topBest.push_back(best);
            Rcpp::Rcout << pos << " score: " << best.object << " P:" << best.P << " N:" << best.N  << " hypothesis: ";


            //best.node2Specified->onto_ref->printSemanticPattern(&(best.bitset[best.intNode2Specified]));

            for(int io = 0; io < refOntologies->size(); ++io)
            {
                if(!best.bitset[io].none())
                {
                    Rcpp::Rcout << (*refOntologies)[io]->getName() << "_";
                    (*refOntologies)[io]->printSemanticPattern(&(best.bitset[io]));
                }
            }

            Rcpp::Rcout << std::endl;
            pos++;
        }

        for(int ii = 0; ii < topBest.size(); ++ii)
        {
            max_heap_TOTAL.push(topBest[ii]);
        }

        OPEN_heap.push(max_heap_TOTAL.top());
        Rcpp::Rcout << "cycle END " << icycles << std::endl << std::endl << std::endl;
    }

    //next steps
    //prepare CLOSED
//    for(int i = 0; i < refOntologies->size(); ++i)
//    {
//        CLOSED[i].clear();
//    }
//    while(!OPEN_heap.empty())
//    {
//        OPEN_heap.pop();
//    }
    //generate new candidates

    //conjunction_max bb = max_heap.top();


    /*
     * pro k nejlepsich termu vygeneruj nova pravidla, tj. pridej ke kazdemu roots
     * */
//    std::vector<std::priority_queue<conjunction_max2> > toSearch = getElongatedRules(&max_heap);

    /*
     * pro kazde takove pravidlo najdi k nejlepsich
     * skore k nejlepsich porovnej se skorem predchazejiciho a vyber jen ty lepsi, horsi zahod
     * */


    /*
     * opakuj dokud nedosahnes pozadovane maximalni delky pravidla
     */



/*

    std::cout << "Stats: best " << QUEUE_LIMIT << " candidates" << std::endl;
    int pos = 1;
    while(pos != 11)//!max_heap.empty())
    {
        conjunction_max best = max_heap.top();
        max_heap.pop();
        std::cout << pos << " score: " << best.object << " hypothesis: ";
        best.node2Specified->onto_ref->printSemanticPattern(&(best.bitset[best.intNode2Specified]));
        pos++;
    }
    */

    return 0;
}

int BeamTopDown::runTest()
{
    std::priority_queue<conjunction_max> max_heap;
    std::priority_queue<conjunction_max> max_heap_TOTAL;
    std::priority_queue<conjunction_min> min_heap;
    std::priority_queue<conjunction_max> OPEN_heap;

    int BEAM_LIMIT = 300;

    //CLOSED for each ontology
    std::vector<boost::dynamic_bitset<> > CLOSED;
    for(int i = 0; i < refOntologies->size(); ++i)
    {
        boost::dynamic_bitset<> ontoBit((*refOntologies)[i]->getOntologyParser()->getTermsNumber(), 0);
        CLOSED.push_back(ontoBit);
    }
    conjunction_max bestConjunction;
    //first step
    initPriorityQueue(&OPEN_heap, &max_heap, &CLOSED);
    conjunction_max tmpBest = OPEN_heap.top();
    OPEN_heap.pop();

    max_heap_TOTAL = OPEN_heap;

    while(!OPEN_heap.empty())
    {
        OPEN_heap.pop();
    }
    OPEN_heap.push(tmpBest);

   // initPriorityQueue(&OPEN_heap, &max_heap, &CLOSED);
    for(int i = 0; i < ITERATE_LIMIT; ++i)
    {

            if(OPEN_heap.empty())
                    break;

            generateNewConjunctions_COMPLETE(&OPEN_heap, &max_heap, &CLOSED);
            //if(!QUEUE_UNLIMITED)
            //        checkQueueMaxLimit(&max_heap, &min_heap, 100000);
    }

    Rcpp::Rcout << "Stats: best " << QUEUE_LIMIT << " candidates" << std::endl;

    conjunction_max actbest = max_heap.top();
    max_heap.pop();

    Rcpp::Rcout << "Best of ALL: ";
    Rcpp::Rcout << "onto: " << actbest.node2Specified->onto_ref->getName() << " score: " << actbest.object << " P:" << actbest.P << " N:" << actbest.N  << " hypothesis: ";
    actbest.node2Specified->onto_ref->printSemanticPattern(&(actbest.bitset[actbest.intNode2Specified]));
    Rcpp::Rcout << std::endl;

    return 0;
}

bool swap(std::priority_queue<conjunction_max> *max_heap, std::priority_queue<conjunction_min> *min_heap)
{
        while(!max_heap->empty())
        {
                conjunction_max tmp = max_heap->top();
                conjunction_min tmp_min;

                tmp_min.bitset = tmp.bitset;
                tmp_min.intNode2Specified = tmp.intNode2Specified;
                tmp_min.object = tmp.object;
                tmp_min.ontoID = tmp.ontoID;
                tmp_min.node2Specified = tmp.node2Specified;
                tmp_min.discovered = tmp.discovered;
                tmp_min.P = tmp.P;
                tmp_min.N = tmp.N;
                tmp_min.bestPathScore = tmp.bestPathScore;

                min_heap->push(tmp_min);
                max_heap->pop();
        }
}

bool swap(std::priority_queue<conjunction_min> *min_heap, std::priority_queue<conjunction_max> *max_heap)
{
        while(!min_heap->empty())
        {
                conjunction_min tmp = min_heap->top();
                conjunction_max tmp_max;

                tmp_max.bitset = tmp.bitset;
                tmp_max.intNode2Specified = tmp.intNode2Specified;
                tmp_max.object = tmp.object;
                tmp_max.ontoID = tmp.ontoID;
                tmp_max.node2Specified = tmp.node2Specified;
                tmp_max.discovered = tmp.discovered;
                tmp_max.P = tmp.P;
                tmp_max.N = tmp.N;
                tmp_max.bestPathScore = tmp.bestPathScore;

                max_heap->push(tmp_max);
                min_heap->pop();
        }
}

void BeamTopDown::checkQueueMaxLimit(std::priority_queue<conjunction_max> *max_heap, std::priority_queue<conjunction_min> *min_heap, int limit)
{
        int heap_size = max_heap->size();
        if(heap_size > limit)
        {
                int toDelete = heap_size - limit;
                swap(max_heap, min_heap);
                for(int i = 0; i <= toDelete; ++i)
                {
                        if(!min_heap->empty())
                        {
                            //conjunction_min toDel = min_heap->top();
                            //delete toDel.bitset;
                            min_heap->pop();
                        }
                }
                swap(min_heap, max_heap);
        }
}

void BeamTopDown::initPriorityQueue(std::priority_queue<conjunction_max> *OPEN_heap, std::priority_queue<conjunction_max> *max_heap, std::vector<boost::dynamic_bitset<> > *CLOSED)
{
    for(int iOnto = 0; iOnto < this->refOntologies->size(); ++iOnto)
    {
        std::vector<Node*> *roots = (*refOntologies)[iOnto]->getOntologyParser()->getRoots();
        //from bottom up to TOP DOWN
        for(int i = 0; i < roots->size(); ++i)
        {

            (*roots)[i]= (*roots)[i]->onto_ref->getOntologyParser()->getNodeTopDown(&((*roots)[i]->id));
        }

        for(int iroot = 0; iroot < roots->size(); ++iroot)
        {
            //initialization
            conjunction_max newItem;
            newItem.node2Specified = (*roots)[iroot];   //set current node!
            newItem.bitset = std::vector<boost::dynamic_bitset<> >(this->refOntologies->size());
            newItem.ontoID = std::vector<int>(this->refOntologies->size(),0);
            for(int ihypothesis = 0; ihypothesis < this->refOntologies->size(); ++ihypothesis)
            {
                newItem.ontoID[ihypothesis] = ihypothesis;
                boost::dynamic_bitset<> bitset((*refOntologies)[ihypothesis]->getOntologyParser()->getTermsNumber(),0);
                newItem.bitset[ihypothesis] = bitset;
            }
            newItem.discovered = newItem.bitset;
            newItem.bitset[iOnto][(*roots)[iroot]->idNum]  = 1; //set root in the current ontology
            newItem.intNode2Specified = iOnto;
            newItem.object = this->evaluateF1Score(&newItem);//newItem.bitset[ontoID], &ontoID);
            newItem.bestPathScore = newItem.object;

            Rcpp::Rcout << "root obj: " << newItem.object  << " ONTO: " << (*roots)[iroot]->onto_ref->getName() << " P: " << newItem.P << " N: " << newItem.N << std::endl;
            (*CLOSED)[iOnto][(*roots)[iroot]->idNum] = 1;//set closed

            max_heap->push(newItem);
            OPEN_heap->push(newItem);
        }
    }
}


void BeamTopDown::initPriorityQueue2(conjunction_max *toInit, std::priority_queue<conjunction_max> *OPEN_heap, std::priority_queue<conjunction_max> *max_heap, std::vector<boost::dynamic_bitset<> > *CLOSED)
{
    for(int iOnto = 0; iOnto < this->refOntologies->size(); ++iOnto)
    {
        std::vector<Node*> *roots = (*refOntologies)[iOnto]->getOntologyParser()->getRoots();
        //from bottom up to TOP DOWN
        for(int i = 0; i < roots->size(); ++i)
        {
            (*roots)[i]= (*roots)[i]->onto_ref->getOntologyParser()->getNodeTopDown(&((*roots)[i]->id));
        }
    }
}


double BeamTopDown::evaluateF1Score(boost::dynamic_bitset<> *bitset, int *ontoID)
{
    int tp = 0;
    int fp = 0;
    int fn = 0;
    for(int iposExample = 0; iposExample < (*this->POSexamples).size(); ++iposExample)
    {
        boost::dynamic_bitset<> *bitsetExample = (*this->POSexamples)[iposExample][*ontoID];//(*this->POSexamples)[*ontoID][iposExample];
        if((*bitset & *bitsetExample).any())
            tp++;
        else
            fn++;
    }

    for(int inegExample = 0; inegExample < (*this->NEGexamples).size(); ++inegExample)
    {
        boost::dynamic_bitset<> *bitsetExample = (*this->NEGexamples)[inegExample][*ontoID];
        if((*bitset & *bitsetExample).any())
            fp++;
    }

    return (2*tp)/(double)(2*tp + fp + fn);
}

double BeamTopDown::evaluateF1Score(conjunction_max *item)
{
    int tp = 0;
    int fp = 0;
    int fn = 0;
    for(int iposExample = 0; iposExample < (*this->POSexamples).size(); ++iposExample)
    {
        bool isCovered = false;
        for(int ibitset = 0; ibitset < item->bitset.size(); ++ibitset)
        {
            //at least one time evaluated
            if((*item).bitset[ibitset].any())
            {
                //example is covered
                boost::dynamic_bitset<> *bitsetExample = (*this->POSexamples)[iposExample][ibitset];
                if(((*item).bitset[ibitset] & *bitsetExample) == (*item).bitset[ibitset])
                {
                    isCovered = true;
                }
                else
                {
                    isCovered = false;
                    break;
                }
            }
        }
        if(isCovered)
            tp++;
        else
            fn++;
    }

    if(tp == 0)
    {
        item->P = 0;
        item->N = 0;
        return 0;
    }

    for(int inegExample = 0; inegExample < (*this->NEGexamples).size(); ++inegExample)
    {
        bool isCovered = false;
        for(int ibitset = 0; ibitset < item->bitset.size(); ++ibitset)
        {
            //at least one time evaluated
            if((*item).bitset[ibitset].any())
            {
                boost::dynamic_bitset<> *bitsetExample = (*this->NEGexamples)[inegExample][ibitset];
                //example is covered
                if(((*item).bitset[ibitset] & *bitsetExample) == (*item).bitset[ibitset])
                {
                    isCovered = true;
                }
                else
                {
                    isCovered = false;
                    break;
                }
            }
        }

        if(isCovered)
            fp++;
    }


    item->P = tp;
    item->N = fp;

    return (2*tp)/(double)(2*tp + fp + fn);
}

bool BeamTopDown::generateNewConjunctions(std::priority_queue<conjunction_max> *OPEN_heap, std::priority_queue<conjunction_max> *max_heap, std::vector<boost::dynamic_bitset<> > *CLOSED)//, vector<boost::dynamic_bitset<> > *featureBitSet, boost::unordered_map<string, int> *CLOSED, boost::dynamic_bitset<> *classMask, std::vector<conjunction_max> *baseTerms)
{
    conjunction_max best = OPEN_heap->top();
    OPEN_heap->pop();
    int debug = 0;    

    Node *toSpecified = best.node2Specified;

    //generate ancestors
    for(std::vector<edge>::iterator iedge = toSpecified->relationship.begin(); iedge != toSpecified->relationship.end(); ++iedge)
    {

        if(iedge->end != NULL)
        {
            //CLOSED test
            if(best.discovered[best.intNode2Specified][iedge->end->idNum] == 1)
            {
                conjunction_max sibling = best;
                sibling.node2Specified = iedge->end;
                OPEN_heap->push(sibling);
            }
            //not in CLOSED list or has not been in a hypothesis list
            else if(((*CLOSED)[best.intNode2Specified][iedge->end->idNum]) == 0 && best.bitset[best.intNode2Specified][iedge->end->idNum] == 0)
            {
                (*CLOSED)[best.intNode2Specified][iedge->end->idNum] = 1; //set CLOSED
                (*CLOSED)[best.intNode2Specified][iedge->start->idNum] = 1; //set CLOSED


                conjunction_max sibling;// = best;
                sibling.bitset = best.bitset;
                sibling.intNode2Specified = best.intNode2Specified;
                sibling.node2Specified = best.node2Specified;
                sibling.discovered = best.discovered;
                sibling.object = 0;
                sibling.ontoID = best.ontoID;

                sibling.bitset[sibling.intNode2Specified][iedge->start->idNum] = 0;
                sibling.discovered[sibling.intNode2Specified][iedge->start->idNum] = 1; //already discovered, will be more general

                if(best.node2Specified->idNum != iedge->start->idNum)
                    Rcpp::Rcout << "edge: " << iedge->start->idNum << " best: " << best.node2Specified->idNum << std::endl;

                sibling.bitset[sibling.intNode2Specified][iedge->end->idNum] = 1; //set current
                sibling.node2Specified = iedge->end;

                double actScore = this->evaluateF1Score(&sibling);
                sibling.object = actScore;

                if(actScore > (best.bestPathScore - EPS_DOUBLE))
                    sibling.bestPathScore = actScore;
                else
                    sibling.bestPathScore = best.bestPathScore;

                if(debug)
                {
                    Rcpp::Rcout << "score: " << actScore << " id: " << iedge->end->idNum << " onto id|sib|best: " << sibling.intNode2Specified << best.intNode2Specified <<  " best score:" << best.bestPathScore << " siblingbest: " << sibling.bestPathScore << " parent: " << iedge->start->idNum << " ";
                    sibling.node2Specified->onto_ref->printSemanticPattern(&(sibling.bitset[sibling.intNode2Specified]));
                }

                double potentialBestScore = (2*sibling.P) / (double)(2*sibling.P + (this->POSexamples->size() - sibling.P));
                if(debug)
                    Rcpp::Rcout << " potential: "  << potentialBestScore << " P: " << sibling.P << " N: " << sibling.N;
                if(potentialBestScore > (best.bestPathScore - EPS_DOUBLE))// || actScore > (best.object - EPS_DOUBLE))
                {
                    OPEN_heap->push(sibling);
                    max_heap->push(sibling);
                    if(debug)
                        Rcpp::Rcout << " OK!!";
                }
                bool a = potentialBestScore > (best.bestPathScore - EPS_DOUBLE);
                bool b = actScore > (best.object - EPS_DOUBLE);

                if(debug)
                {
                    Rcpp::Rcout << " - pot cond: " << a << " act cond: " << b;
                    Rcpp::Rcout << std::endl;
                }
            }
        }
    }
}

bool BeamTopDown::generateNewConjunctions_COMPLETE(std::priority_queue<conjunction_max> *OPEN_heap, std::priority_queue<conjunction_max> *max_heap, std::vector<boost::dynamic_bitset<> > *CLOSED)//, vector<boost::dynamic_bitset<> > *featureBitSet, boost::unordered_map<string, int> *CLOSED, boost::dynamic_bitset<> *classMask, std::vector<conjunction_max> *baseTerms)
{
    conjunction_max best = OPEN_heap->top();
    OPEN_heap->pop();
    Node *toSpecified = best.node2Specified;

    //generate ancestors
    for(std::vector<edge>::iterator iedge = toSpecified->relationship.begin(); iedge != toSpecified->relationship.end(); ++iedge)
    {

        if(iedge->end != NULL)
        {
            //CLOSED test
            //not in CLOSED list or has not been in a hypothesis list
            if(!((*CLOSED)[best.intNode2Specified].test(iedge->end->idNum)))
            {
                (*CLOSED)[best.intNode2Specified][iedge->end->idNum] = 1; //set CLOSED


                conjunction_max sibling;// = best;
                sibling.bitset = best.bitset;
                sibling.intNode2Specified = best.intNode2Specified;
                sibling.node2Specified = best.node2Specified;
                sibling.discovered = best.discovered;
                sibling.object = 0;
                sibling.ontoID = best.ontoID;

                sibling.bitset[sibling.intNode2Specified][iedge->start->idNum] = 0;
                sibling.discovered[sibling.intNode2Specified][iedge->start->idNum] = 1; //already discovered, will be more general

                if(best.node2Specified->idNum != iedge->start->idNum)
                    Rcpp::Rcout << "edge: " << iedge->start->idNum << " best: " << best.node2Specified->idNum << std::endl;


                sibling.bitset[sibling.intNode2Specified][iedge->end->idNum] = 1; //set current
                sibling.node2Specified = iedge->end;
                double actScore = this->evaluateF1Score(&sibling);
                sibling.object = actScore;

                Rcpp::Rcout << "score: " << actScore << " id: " << iedge->end->idNum << " hypothesis: ";
                sibling.node2Specified->onto_ref->printSemanticPattern(&(sibling.bitset[sibling.intNode2Specified]));
                Rcpp::Rcout << std::endl;

                            OPEN_heap->push(sibling);
                            max_heap->push(sibling);

            }
        }
    }

}

bool BeamTopDown::isBetter(conjunction_max *bestConjunction, conjunction_max *max_heap_top)
{
        if(bestConjunction->object < max_heap_top->object)
                return true;
        return false;
}

