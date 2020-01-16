#include "example.h"
/*
Example::Example(int kbics, arma::mat *armaData, std::vector<paretoSet> *myparetoSet, std::vector<Ontology *> *refOntologies)
{
    this->kbics = kbics;
    this->armaData = armaData;
    this->myparetoSet = myparetoSet;
    this->refOntologies = refOntologies;
}
*/
Example::Example(int kbics, arma::mat *armaData, std::vector<Ontology *> *refOntologies, int verbose)
{
    this->kbics = kbics;
    this->armaData = armaData;
    this->refOntologies = refOntologies;
    this->verbose = verbose;
}
//for testing
Example::Example(arma::mat *testarmaData, boost::unordered_map<std::string, std::vector<std::vector<std::string> > > *rowTestOntologyDesc, boost::unordered_map<std::string, std::vector<std::vector<std::string> > > *colTestOntologyDesc, std::vector<Ontology *> *refOntologies, int verbose)
{
    this->armaData = testarmaData;
    this->refOntologies = refOntologies;
    this->rowTestOntologyDesc = rowTestOntologyDesc;
    this->colTestOntologyDesc = colTestOntologyDesc;
    this->verbose = verbose;
}


void Example::findPositiveExamples(double threshold)
{
    /*
    arma::mat sumMatOverAll(armaData->n_rows, armaData->n_cols, arma::fill::zeros);
    //pareto bicluster -> example->ontology->bitset
    std::vector< std::vector< std::vector<boost::dynamic_bitset<>* > > > POSexamples;//tPos.n_cols);

    for(int ibic = 0; ibic < this->kbics; ++ibic)
    {
        //semanticBicluster tmpBic;
        arma::mat sumMat(armaData->n_rows, armaData->n_cols, arma::fill::zeros);
        int irow = 0;
        for(int iparetoBic = ibic; iparetoBic < bicAssign->size(); iparetoBic += this->kbics)
        {
            arma::mat tmp = (*myparetoSet)[irow].chromozome[2*ibic] * arma::trans((*myparetoSet)[irow].chromozome[2*ibic+1]);
            sumMat += tmp;
            sumMatOverAll += tmp;
            ++irow;
        }

        //find positive examples in bicluster
        arma::uvec indices = arma::find(sumMat > (threshold * arma::max(arma::max(sumMat))));
        arma::umat tPos = arma::ind2sub(arma::size(sumMat), indices);

        std::cout << "Positive examples: " << tPos.n_cols << std::endl;


        //POSITIVE EXAMPLE
        std::vector< std::vector<boost::dynamic_bitset<>* > > actPositive(tPos.n_cols);
        this->buildTrainingExamples(&actPositive, &tPos, refOntologies, &indexPosExamples);
        POSexamples.push_back(actPositive);
    }
    this->sumMatOverAll = sumMatOverAll;
    this->POSexamples = POSexamples;
    */
}

void Example::findNegativeExamples(double threshold)
{
    /*
    arma::uvec indicesNeg = arma::find(this->sumMatOverAll < (threshold * arma::max(arma::max(this->sumMatOverAll))));
    arma::umat tNeg = arma::ind2sub(arma::size(this->sumMatOverAll), indicesNeg);

    std::vector< std::vector<boost::dynamic_bitset<>* > > NEGexamples(tNeg.n_cols);

    this->buildTrainingExamples(&NEGexamples, &tNeg, refOntologies);
    this->NEGexamples = NEGexamples;

    std::cout << "Negative examples: " << tNeg.n_cols << std::endl;
    */
}

//sem1r version !!!!!!!
void Example::findPositiveExamples()
{
    //arma::mat sumMatOverAll(armaData->n_rows, armaData->n_cols, arma::fill::zeros);
    //example->ontology->bitset
    std::vector< std::vector<boost::dynamic_bitset<>* > > POSexamples;//tPos.n_cols);

    //find positive examples in bicluster
    arma::uvec indices = arma::find(*(this->armaData) == 1);
    arma::umat tPos = arma::ind2sub(arma::size(*(this->armaData)), indices);

    if(this->verbose)
        Rcpp::Rcout << "Positive examples: " << tPos.n_cols << std::endl;


    //POSITIVE EXAMPLE
    std::vector< std::vector<boost::dynamic_bitset<>* > > actPositive(tPos.n_cols);
    this->buildTrainingExamples(&actPositive, &tPos, refOntologies, &indexPosExamples);
    //POSexamples.push_back(actPositive);
    POSexamples = actPositive;

    //this->sumMatOverAll = sumMatOverAll;
    this->POSexamples = POSexamples;
}

void Example::findNegativeExamples()
{
    arma::uvec indicesNeg = arma::find(*(this->armaData) == 0);//this->sumMatOverAll < (threshold * arma::max(arma::max(this->sumMatOverAll))));
    arma::umat tNeg = arma::ind2sub(arma::size(*(this->armaData)), indicesNeg);

    std::vector< std::vector<boost::dynamic_bitset<>* > > NEGexamples(tNeg.n_cols);

    this->buildTrainingExamples(&NEGexamples, &tNeg, refOntologies);
    this->NEGexamples = NEGexamples;

    if(this->verbose)
        Rcpp::Rcout << "Negative examples: " << tNeg.n_cols << std::endl;
}
//sem1r version END !!!!!!!


Example::~Example()
{

}

void Example::buildTrainingExamples(std::vector< std::vector<boost::dynamic_bitset<>* > > *examples, arma::umat *indexMatrix, std::vector<Ontology *> *refOntologies)
{

    //iterate over all potencial examples
    for(int iElem = 0; iElem < indexMatrix->n_cols; ++iElem)
    {
        int row = indexMatrix->col(iElem)[0]; //row index
        int col = indexMatrix->col(iElem)[1];  //col index

        for(int iOnto = 0; iOnto < refOntologies->size(); ++iOnto)
        {
            if((*refOntologies)[iOnto]->getOntoType() == ROW_ONTOLOGY)
                (*examples)[iElem].push_back((*refOntologies)[iOnto]->getSemanticPattern(row));
            else
                (*examples)[iElem].push_back((*refOntologies)[iOnto]->getSemanticPattern(col));
        }
    }
}

void Example::buildTrainingExamples(std::vector< std::vector<boost::dynamic_bitset<>* > > *examples, arma::umat *indexMatrix, std::vector<Ontology *> *refOntologies, std::vector<boost::unordered_map<int, bool> > *indexPosExamples)
{
    //prepare hash of indexes for positive examples
    for(int iOnto = 0; iOnto < refOntologies->size(); ++iOnto)
    {
        boost::unordered_map<int, bool> tmpHash;
        indexPosExamples->push_back(tmpHash);
    }
    //iterate over all potencial examples
    for(int iElem = 0; iElem < indexMatrix->n_cols; ++iElem)
    {
        int row = indexMatrix->col(iElem)[0]; //row index
        int col = indexMatrix->col(iElem)[1];  //col index

        for(int iOnto = 0; iOnto < refOntologies->size(); ++iOnto)
        {
            if((*refOntologies)[iOnto]->getOntoType() == ROW_ONTOLOGY)
            {
                (*examples)[iElem].push_back((*refOntologies)[iOnto]->getSemanticPattern(row));
                (*indexPosExamples)[iOnto][row] = 1;   //row index for positive examples
            }
            else
            {
                (*examples)[iElem].push_back((*refOntologies)[iOnto]->getSemanticPattern(col));
                (*indexPosExamples)[iOnto][col] = 1; //col index for positive examples
            }
        }
    }
}

boost::dynamic_bitset<> Example::coveredExamplesByRule(newComplexStat *rule)
{
    boost::dynamic_bitset<> covered(this->POSexamples.size() + this->NEGexamples.size());
    covered.set();  //set all bit to 1
    //iterate over positive example
    //std::cout << "covered example posotive..." << std::endl;
    for(int ipos = 0; ipos < this->POSexamples.size(); ++ipos)
    {
        //iterate over onto
        for(int ionto = 0; ionto < this->POSexamples[ipos].size(); ++ionto)
        {
            for(int irule = 0; irule < rule->rules.size(); ++irule)
            {
                //this->POSexamples[ipos][ionto] - bitset of terms for ref ontology
                //actrule.rules[iterm]->termRefBottomUp->idNum;

    //            std::cout << "ionto: " << ionto << " ipos: " << ipos << " irule" << irule << std::endl;
    //            std::cout << "name:" << rule->rules[irule]->termRefBottomUp->onto_ref->getName() << std::endl;
 //               if(rule->rules[irule]->noderef != NULL)
                 //   std::cout << "name:" << rule->rules[irule]->noderef->name << std::endl;

                if(rule->rules[irule].noderef->onto_ref->getName() == (*this->refOntologies)[ionto]->getName())
                {
                                //print test rule
                                        //std::cout << "rule: " <<rule->rules[irule].noderef->idNum << " ones: " << this->POSexamples[ipos][ionto]->count() << " size: " << this->POSexamples[ipos][ionto]->size() << std::endl;

                    //if term is positive
                    if(rule->rules[irule].posFeature)
                    {
                        if(this->POSexamples[ipos][ionto]->test(rule->rules[irule].noderef->idNum) == 0)
                        {
                            covered[ipos] = 0;
                        }
                    }
                    else    //negative
                    {
                        if(this->POSexamples[ipos][ionto]->test(rule->rules[irule].noderef->idNum) == 1)
                        {
                            covered[ipos] = 0;
                        }
                    }

                }
                //rule->rules[irule].termRefBottomUp->onto_ref->getSemanticPatternTest(ipos);

            }

        }
    }
    //negative
    int inegIndex = this->POSexamples.size();
    //std::cout << "positive DONE" << std::endl << "covered example negative..." << std::endl;
    for(int ineg = 0; ineg < this->NEGexamples.size(); ++ineg)
    {
        //iterate over onto
        for(int ionto = 0; ionto < this->NEGexamples[ineg].size(); ++ionto)
        {
            for(int irule = 0; irule < rule->rules.size(); ++irule)
            {
                //this->POSexamples[ipos][ionto] - bitset of terms for ref ontology
                //actrule.rules[iterm]->termRefBottomUp->idNum;
                if(rule->rules[irule].noderef->onto_ref->getName() == (*this->refOntologies)[ionto]->getName())
                {
                    if(rule->rules[irule].posFeature)
                    {
                        if(this->NEGexamples[ineg][ionto]->test(rule->rules[irule].noderef->idNum) == 0)
                        {
                            covered[inegIndex] = 0;
                        }
                    }
                    else
                    {
                        if(this->NEGexamples[ineg][ionto]->test(rule->rules[irule].noderef->idNum) == 1)
                        {
                            covered[inegIndex] = 0;
                        }
                    }
                }
                //rule->rules[irule]->termRefBottomUp->onto_ref->getSemanticPatternTest(ipos)

            }

        }
        ++inegIndex;
    }
    return covered;
}



void Example::buildTestingExamples(double positiveTH, double negativeTH)
{
    //find positive examples in bicluster
    arma::uvec indices = arma::find(*(this->armaData) == 1);
    arma::umat tPos = arma::ind2sub(arma::size(*(this->armaData)), indices);

    if(this->verbose)
        Rcpp::Rcout << "Positive examples: " << tPos.n_cols << std::endl;


    //POSITIVE EXAMPLE
    //std::vector< std::vector<boost::dynamic_bitset<>* > > actPositive(tPos.n_cols);
    //this->buildTrainingExamples(&actPositive, &tPos, refOntologies, &indexPosExamples);

    //negative example
    arma::uvec indicesNeg = arma::find(*(this->armaData) == 0);//this->sumMatOverAll < (threshold * arma::max(arma::max(this->sumMatOverAll))));
    arma::umat tNeg = arma::ind2sub(arma::size(*(this->armaData)), indicesNeg);

    if(this->verbose)
        Rcpp::Rcout << "Negative examples: " << tNeg.n_cols << std::endl;

    std::vector< std::vector<boost::dynamic_bitset<>* > > POSexamples(tPos.n_cols);
    std::vector< std::vector<boost::dynamic_bitset<>* > > NEGexamples(tNeg.n_cols);

    //this->buildTrainingExamples(&NEGexamples, &tNeg, refOntologies);


    //prepare hash of indexes for positive examples
    //for(int iOnto = 0; iOnto < refOntologies->size(); ++iOnto)
    //{
    //    boost::unordered_map<int, bool> tmpHash;
    //    indexPosExamples->push_back(tmpHash);
    //}
    //iterate over all potencial examples
    std::vector< std::vector<boost::dynamic_bitset<>* > > examples(tPos.n_cols + tNeg.n_cols);

    this->classBitMask.resize(tPos.n_cols + tNeg.n_cols, 0);

//    std::cout << "1" << std::endl;
//    std::cout << "positive..." ;
    boost::dynamic_bitset<> myclass(tPos.n_cols + tNeg.n_cols, 0);
    for(int iElem = 0; iElem < tPos.n_cols; ++iElem)
    {
        int row = tPos.col(iElem)[0]; //row index
        int col = tPos.col(iElem)[1];  //col index
        //std::cout << "col: " << col << " row: " << row << std::endl;
//std::cout << "2" << std::endl;
        //std::cout << "onto start.." << std::endl;
        for(int iOnto = 0; iOnto < refOntologies->size(); ++iOnto)
        {
            //std::cout << "3" << std::endl;
//            std::cout << "i: " << iOnto << std::endl;
            if((*refOntologies)[iOnto]->getOntoType() == ROW_ONTOLOGY)
            {
                //std::cout << "4" << std::endl;
                POSexamples[iElem].push_back((*refOntologies)[iOnto]->getSemanticPatternTest(row));
                examples[iElem].push_back((*refOntologies)[iOnto]->getSemanticPatternTest(row));
                //(*indexPosExamples)[iOnto][row] = 1;   //row index for positive examples
                //std::cout << "5" << std::endl;
            }
            else
            {
                //std::cout << "6" << std::endl;
                POSexamples[iElem].push_back((*refOntologies)[iOnto]->getSemanticPatternTest(col));
                examples[iElem].push_back((*refOntologies)[iOnto]->getSemanticPatternTest(col));
                //(*indexPosExamples)[iOnto][col] = 1; //col index for positive examples
                //std::cout << "7" << std::endl;
            }
        }
        //std::cout << "8" << std::endl;
        //std::cout << "onto end" << std::endl;
        myclass[iElem] = 1;
        this->classBitMask[iElem] = 1;
    }
//std::cout << "DONE" << std::endl << "negative...";
    for(int iElem = 0; iElem < tNeg.n_cols; ++iElem)
    {
        int row = tNeg.col(iElem)[0]; //row index
        int col = tNeg.col(iElem)[1];  //col index

        for(int iOnto = 0; iOnto < refOntologies->size(); ++iOnto)
        {
            if((*refOntologies)[iOnto]->getOntoType() == ROW_ONTOLOGY)
            {
                NEGexamples[iElem].push_back((*refOntologies)[iOnto]->getSemanticPatternTest(row));
                examples[iElem].push_back((*refOntologies)[iOnto]->getSemanticPatternTest(row));
                //(*indexPosExamples)[iOnto][row] = 1;   //row index for positive examples
            }
            else
            {
                NEGexamples[iElem].push_back((*refOntologies)[iOnto]->getSemanticPatternTest(col));
                examples[iElem].push_back((*refOntologies)[iOnto]->getSemanticPatternTest(col));
                //(*indexPosExamples)[iOnto][col] = 1; //col index for positive examples
            }
        }
    }
//std::cout << "DONE" << std::endl;
    this->POSexamples = POSexamples;
    this->NEGexamples = NEGexamples;
}


Rcpp::DataFrame Example::getPropositionalTable(int bic)
{
    Rcpp::DataFrame propositionTable;
    int featureNumber = 0;
    std::vector<std::string> featureName;
    for(int ife = 0; ife < refOntologies->size(); ++ife)
    {
        featureNumber += (*refOntologies)[ife]->getOntologyParser()->getTermsNumber();
        for(int ibit = 0; ibit < (*refOntologies)[ife]->getOntologyParser()->getTermsNumber(); ++ibit)
        {
             featureName.push_back((*refOntologies)[ife]->getOntologyParser()->getNodesByPosition(ibit)->name);
        }
    }
    //->feature -> example
    std::vector<Rcpp::LogicalVector> featureVector(featureNumber);
    //featureVector.resize(featureNumber);
    Rcpp::LogicalVector featureClass(POSexamples.size() + NEGexamples.size());
    //featureClass.resize(POSexamples[0].size() + NEGexamples.size());

    if(this->verbose)
        Rcpp::Rcout << "feature number: " << featureNumber << std::endl;
    //initialize vector
    for(int iinit = 0; iinit < featureVector.size(); ++iinit)
    {
        Rcpp::LogicalVector tmp(POSexamples.size() + NEGexamples.size());
        featureVector[iinit] = tmp;
    }

    //bicluster 0 - positive examples
    for(int iexpos = 0; iexpos < POSexamples.size(); ++iexpos)
    {
        int ivect = 0;
        for(int ionto = 0; ionto < POSexamples[iexpos].size(); ++ionto)
        {
            for(int ibit = 0; ibit < POSexamples[iexpos][ionto]->size(); ++ibit)
            {
                if(POSexamples[iexpos][ionto]->test(ibit))
                    featureVector[ibit+ivect][iexpos] = true;
                else
                    featureVector[ibit+ivect][iexpos] = false;
            }
            ivect += POSexamples[iexpos][ionto]->size();
        }
        featureClass[iexpos] = true;
    }

    //bicluster 0 - negative examples
    for(int iexpos = 0; iexpos < NEGexamples.size(); ++iexpos)
    {
        int ivect = 0;
        for(int ionto = 0; ionto < NEGexamples[iexpos].size(); ++ionto)
        {
            for(int ibit = 0; ibit < NEGexamples[iexpos][ionto]->size(); ++ibit)
            {
                if(NEGexamples[iexpos][ionto]->test(ibit))
                    featureVector[ibit+ivect][iexpos+POSexamples.size()] = true;
                else
                    featureVector[ibit+ivect][iexpos+POSexamples.size()] = false;
            }
            ivect += NEGexamples[iexpos][ionto]->size();
        }
        featureClass[iexpos + POSexamples.size()] = false;
    }


    for(int it = 0; it < featureVector.size(); ++it)
    {
        //at least one covered
        bool isCovered = false;
        for(int ifv = 0; ifv < featureVector[it].size(); ++ifv)
        {
            if(featureVector[it][ifv] == true)
            {
                isCovered = true;
                break;
            }
        }
        if(isCovered)
        {
                propositionTable.push_back(featureVector[it], featureName[it]);
        }
    }
    propositionTable.push_back(featureClass, "Class");
    return propositionTable;
}

Rcpp::DataFrame Example::getEnrichTable(int bic, std::vector<Node *> *enrichNodes)
{
    Rcpp::DataFrame enrichPropositionTable;
    std::vector<boost::dynamic_bitset<> > enrichItems;
    std::vector<boost::unordered_map<unsigned int, bool> > enrichIDnums;

    int featureNumber = 0;
    for(int ionto = 0; ionto < refOntologies->size(); ++ionto)
    {
        int sizePop = 0;
        int sampleSize = 0;

        int nterm = (*refOntologies)[ionto]->getOntologyParser()->getTermsNumber();
        Rcpp::NumericVector termPosSucc(nterm, 0.0);
        Rcpp::NumericVector termNegSucc(nterm, 0.0);
        Rcpp::NumericVector pVAL(nterm, 0.0);

        Ontology *actOnto = (*refOntologies)[ionto];
        std::vector<boost::dynamic_bitset<> *> patterns = actOnto->getSemanticPatterns();
        int ppos = 0;
        int nneg = 0;

        for(int iexample = 0; iexample < patterns.size(); ++iexample)
        {


            //background examples (population)
            //if(indexPosExamples[bic][ionto].find(iexample) == indexPosExamples[bic][ionto].end())
            if(indexPosExamples[ionto].find(iexample) == indexPosExamples[ionto].end())
            {
                //size_t index = patterns[iexample]->find_first();
                //while(index != boost::dynamic_bitset<>::npos)
                //{
                //    termNegSucc(index) =  termNegSucc(index) + 1;
                //    index = patterns[iexample]->find_next(index);
                //}
                for(int ibit = 0; ibit < patterns[iexample]->size(); ++ibit)
                {
                    if((*patterns[iexample])[ibit] == 1)
                    {
                        termNegSucc(ibit) = termNegSucc(ibit) + 1;
                    }
                }
            }
            else // positive examples (sample)
            {
                //size_t index = patterns[iexample]->find_first();
                //while(index != boost::dynamic_bitset<>::npos)
                //{
                //    termPosSucc(index) =  termPosSucc(index) + 1;
                //    index = patterns[iexample]->find_next(index);
                //}

                for(int ibit = 0; ibit < patterns[iexample]->size(); ++ibit)
                {
                    if((*patterns[iexample])[ibit] == 1)
                    {
                        termPosSucc(ibit) = termPosSucc(ibit) + 1;
                    }
                }
                ++sampleSize;
            }

            if(patterns[iexample]->any())
                ++sizePop;
        }

        //compute variables
        for(int ivec = 0; ivec < termPosSucc.size(); ++ivec)
        {
            int succSample = termPosSucc(ivec);
            int succPop = termNegSucc(ivec);
            int failurePOP = sizePop-succPop;
            pVAL(ivec) = getHyperGeometricScore(&succSample, &succPop, &failurePOP, &sampleSize);
        }
        //multiple hypothesiss correction - bonferroni
        /* In the case of GO::TermFinder, the value used for the Bonferroni correction is the number of nodes
         *  to which the genes of interest are collectively annotated, excluding those nodes which only have
         *  a single annotation in the background distribution, which a priori cannot be significantly enriched.
        */

        //terms that are collectively annotated
        //arma::uvec indeces = find(termPosSucc == sampleSize);

        //boost::dynamic_bitset<> discovered(nterm);
        int bonferroni = 0;
        for(int ivec = 0; ivec < termPosSucc.size(); ++ivec)
        {
            if(termPosSucc(ivec) >= 1 && termNegSucc(ivec) != 1)
                ++bonferroni;
        }

        double alpha = 0.05;
        if(this->verbose)
        {
            Rcpp::Rcout << std::endl;
            Rcpp::Rcout << "alpha: " << alpha << " number of hypothesis: " << bonferroni << " new alpha: " << alpha/(double)bonferroni  << "sample size: " << sampleSize << std::endl;
        }
        boost::dynamic_bitset<> enrichment(nterm);
        boost::unordered_map<unsigned int, bool> idNumHash;
        int nulla = 0;
        for(int ivec = 0; ivec < pVAL.size(); ++ivec)
        {
            if(pVAL(ivec) <= alpha/(double)bonferroni)
            {
                enrichment[ivec] = 1;
                idNumHash[ivec] = true;
                //refOntologies[ionto].getOntologyParser();
                enrichNodes->push_back((*refOntologies)[ionto]->getOntologyParser()->getNodesByPosition(ivec));
                if(this->verbose)
                    Rcpp::Rcout << "ench score: " << pVAL(ivec) << " node: " << (*refOntologies)[ionto]->getOntologyParser()->getNodesByPosition(ivec)->name << " termPos: " << termPosSucc(ivec) << " termNeg: " << termNegSucc(ivec) << std::endl;
                ++featureNumber;
            }

            if(termPosSucc(ivec) == 0)
                nulla++;
        }
        if(this->verbose)
            Rcpp::Rcout << "enrichment for : " << actOnto->getName() << " en count: " << enrichment.count() << " total: " << enrichment.size() << " pocet null:" << nulla << std::endl;

        enrichItems.push_back(enrichment);
        enrichIDnums.push_back(idNumHash);//refOntologies[ionto]->getOntologyParser()->getNodesByPosition(ibit)->idNum);

    }
    std::vector<Rcpp::LogicalVector> featureVector(featureNumber);
    for(int iinit = 0; iinit < featureVector.size(); ++iinit)
    {
        Rcpp::LogicalVector tmp(POSexamples[bic].size() + NEGexamples.size());
        featureVector[iinit] = tmp;
    }


    Rcpp::LogicalVector featureClass(POSexamples[bic].size() + NEGexamples.size());
    std::vector<std::string> featureName(featureNumber);
    this->computeFeatureVectors(bic, &enrichItems, &featureVector, &featureClass, &featureName);

    for(int ifeat = 0; ifeat < featureVector.size(); ++ifeat)
    {
        enrichPropositionTable.push_back(featureVector[ifeat], featureName[ifeat]);
    }
    enrichPropositionTable.push_back(featureClass, "Class");

    return enrichPropositionTable;
}

Rcpp::DataFrame Example::getEnrichTablev2(int bic, std::vector<Node *> *enrichNodes)
{
    Rcpp::DataFrame enrichPropositionTable;
    std::vector<boost::dynamic_bitset<> > enrichItems;
    std::vector<boost::unordered_map<unsigned int, bool> > enrichIDnums;

    int featureNumber = 0;
    for(int ionto = 0; ionto < refOntologies->size(); ++ionto)
    {
        int sizePop = 0;
        int sampleSize = 0;

        int nterm = (*refOntologies)[ionto]->getOntologyParser()->getTermsNumber();
        Rcpp::NumericVector termPosSucc(nterm, 0.0);
        Rcpp::NumericVector termNegSucc(nterm, 0.0);
        Rcpp::NumericVector pVAL(nterm, 0.0);

        Ontology *actOnto = (*refOntologies)[ionto];
        //std::vector<boost::dynamic_bitset<> *> patterns = actOnto->getSemanticPatterns();
        int ppos = 0;
        int nneg = 0;

        //pareto bicluster -> example->ontology->bitset
        this->POSexamples;
        this->NEGexamples;

        //positive examples
        for(int iexample = 0; iexample < this->POSexamples.size(); ++iexample)
        {
            boost::dynamic_bitset<> *patterns = this->POSexamples[iexample][ionto];//actOnto->getSemanticPatterns();


            for(int ibit = 0; ibit < patterns->size(); ++ibit)
            {
                if((*patterns)[ibit] == 1)
                {
                    termPosSucc(ibit) = termPosSucc(ibit) + 1;
                }
            }
            if(patterns->any())
                ++sampleSize;
        }

        //negative examples
        for(int iexample = 0; iexample < this->NEGexamples.size(); ++iexample)
        {
            boost::dynamic_bitset<> *patterns = this->NEGexamples[iexample][ionto];//actOnto->getSemanticPatterns();
            for(int ibit = 0; ibit < patterns->size(); ++ibit)
            {
                if((*patterns)[ibit] == 1)
                {
                    termNegSucc(ibit) = termNegSucc(ibit) + 1;
                }
            }
            if(patterns->any())
                ++sizePop;
        }

        //compute variables
        for(int ivec = 0; ivec < termPosSucc.size(); ++ivec)
        {
            int succSample = termPosSucc(ivec);
            int succPop = termNegSucc(ivec);
            int failurePOP = sizePop-succPop;
            pVAL(ivec) = getHyperGeometricScore(&succSample, &succPop, &failurePOP, &sampleSize);
        }
        //multiple hypothesiss correction - bonferroni
        /* In the case of GO::TermFinder, the value used for the Bonferroni correction is the number of nodes
         *  to which the genes of interest are collectively annotated, excluding those nodes which only have
         *  a single annotation in the background distribution, which a priori cannot be significantly enriched.
        */

        //terms that are collectively annotated
        //arma::uvec indeces = find(termPosSucc == sampleSize);

        //boost::dynamic_bitset<> discovered(nterm);
        int bonferroni = 0;
        for(int ivec = 0; ivec < termPosSucc.size(); ++ivec)
        {
            if(termPosSucc(ivec) >= 1 && termNegSucc(ivec) != 1)
                ++bonferroni;
        }

        double alpha = 0.05;
        if(this->verbose)
        {
            Rcpp::Rcout << std::endl;
            Rcpp::Rcout << "alpha: " << alpha << " number of hypothesis: " << bonferroni << " new alpha: " << alpha/(double)bonferroni  << "sample size: " << sampleSize << std::endl;
        }

        boost::dynamic_bitset<> enrichment(nterm);
        boost::unordered_map<unsigned int, bool> idNumHash;
        int nulla = 0;
        for(int ivec = 0; ivec < pVAL.size(); ++ivec)
        {
            if(pVAL(ivec) <= alpha/(double)bonferroni)
            {
                enrichment[ivec] = 1;
                idNumHash[ivec] = true;
                //refOntologies[ionto].getOntologyParser();
                enrichNodes->push_back((*refOntologies)[ionto]->getOntologyParser()->getNodesByPosition(ivec));
                //std::cout << "ench score: " << pVAL(ivec) << " node: " << (*refOntologies)[ionto]->getOntologyParser()->getNodesByPosition(ivec)->name << " termPos: " << termPosSucc(ivec) << " termNeg: " << termNegSucc(ivec) << std::endl;
                ++featureNumber;
            }

            if(termPosSucc(ivec) == 0)
                nulla++;
        }
        if(this->verbose)
            Rcpp::Rcout << "enrichment for : " << actOnto->getName() << " en count: " << enrichment.count() << " total: " << enrichment.size() << " pocet null:" << nulla << std::endl;

        enrichItems.push_back(enrichment);
        enrichIDnums.push_back(idNumHash);//refOntologies[ionto]->getOntologyParser()->getNodesByPosition(ibit)->idNum);

    }
    std::vector<Rcpp::LogicalVector> featureVector(featureNumber);
    for(int iinit = 0; iinit < featureVector.size(); ++iinit)
    {
        Rcpp::LogicalVector tmp(POSexamples.size() + NEGexamples.size());
        featureVector[iinit] = tmp;
    }


    Rcpp::LogicalVector featureClass(POSexamples.size() + NEGexamples.size());
    std::vector<std::string> featureName(featureNumber);
    this->computeFeatureVectors(bic, &enrichItems, &featureVector, &featureClass, &featureName);

    for(int ifeat = 0; ifeat < featureVector.size(); ++ifeat)
    {
        enrichPropositionTable.push_back(featureVector[ifeat], featureName[ifeat]);
    }
    enrichPropositionTable.push_back(featureClass, "Class");

    return enrichPropositionTable;
}


mydataframe Example::getEnrichTablev2bitset(int bic, std::vector<Node *> *enrichNodes)
{	
    std::vector<boost::dynamic_bitset<> > enrichPropositionTable;
    std::vector<boost::dynamic_bitset<> > enrichItems;
    //std::vector<boost::unordered_map<unsigned int, bool> > enrichIDnums;

    int featureNumber = 0;
    for(int ionto = 0; ionto < refOntologies->size(); ++ionto)
    {
        int sizePop = 0;
        int sampleSize = 0;

        int nterm = (*refOntologies)[ionto]->getOntologyParser()->getTermsNumber();
        std::vector<int> termPosSucc(nterm, 0);
        std::vector<int> termNegSucc(nterm, 0);
        std::vector<double> pVAL(nterm, 0);

        Ontology *actOnto = (*refOntologies)[ionto];
        //std::vector<boost::dynamic_bitset<> *> patterns = actOnto->getSemanticPatterns();
        int ppos = 0;
        int nneg = 0;

        //pareto bicluster -> example->ontology->bitset
        this->POSexamples;
        this->NEGexamples;
        
        boost::dynamic_bitset<> mynewclass(this->POSexamples.size() + this->NEGexamples.size(), 0);

        //positive examples
        for(int iexample = 0; iexample < this->POSexamples.size(); ++iexample)
        {
            boost::dynamic_bitset<> *patterns = this->POSexamples[iexample][ionto];//actOnto->getSemanticPatterns();


            for(int ibit = 0; ibit < patterns->size(); ++ibit)
            {
                if((*patterns)[ibit] == 1)
                {
                    termPosSucc[ibit] = termPosSucc[ibit] + 1;
                }
            }
            if(patterns->any())
                ++sampleSize;
        }

        //negative examples
        for(int iexample = 0; iexample < this->NEGexamples.size(); ++iexample)
        {
            boost::dynamic_bitset<> *patterns = this->NEGexamples[iexample][ionto];//actOnto->getSemanticPatterns();
            for(int ibit = 0; ibit < patterns->size(); ++ibit)
            {
                if((*patterns)[ibit] == 1)
                {
                    termNegSucc[ibit] = termNegSucc[ibit] + 1;
                }
            }
            if(patterns->any())
                ++sizePop;
        }

        //compute variables
        for(int ivec = 0; ivec < termPosSucc.size(); ++ivec)
        {
            int succSample = termPosSucc[ivec];
            int succPop = termNegSucc[ivec];
            int failurePOP = sizePop-succPop;
            pVAL[ivec] = getHyperGeometricScore(&succSample, &succPop, &failurePOP, &sampleSize);
            std::cout << "pval of " << ivec << " is: " << pVAL[ivec] << " succSample: " << succSample << " succPop: " << succPop <<  " failurePOP: " << failurePOP << " sampleSize: " << sampleSize << " sizePop:" << sizePop  << std::endl;
        }
        //multiple hypothesiss correction - bonferroni
        /* In the case of GO::TermFinder, the value used for the Bonferroni correction is the number of nodes
         *  to which the genes of interest are collectively annotated, excluding those nodes which only have
         *  a single annotation in the background distribution, which a priori cannot be significantly enriched.
        */

        //terms that are collectively annotated
        //arma::uvec indeces = find(termPosSucc == sampleSize);

        //boost::dynamic_bitset<> discovered(nterm);
        int bonferroni = 0;
        for(int ivec = 0; ivec < termPosSucc.size(); ++ivec)
        {
            if(termPosSucc[ivec] >= 1 && termNegSucc[ivec] != 1)
                ++bonferroni;
        }

        std::cout << std::endl;
        double alpha = 0.1;
        std::cout << "alpha: " << alpha << " number of hypothesis: " << bonferroni << " new alpha: " << alpha/(double)bonferroni  << "sample size: " << sampleSize << std::endl;

        boost::dynamic_bitset<> enrichment(nterm);
//        boost::unordered_map<unsigned int, bool> idNumHash;
        int nulla = 0;
        for(int ivec = 0; ivec < pVAL.size(); ++ivec)
        {
            if(pVAL[ivec] <= alpha/(double)bonferroni)
            {
                enrichment[ivec] = 1;
 //               idNumHash[ivec] = true;
                //refOntologies[ionto].getOntologyParser();
                enrichNodes->push_back((*refOntologies)[ionto]->getOntologyParser()->getNodesByPosition(ivec));
                //std::cout << "ench score: " << pVAL(ivec) << " node: " << (*refOntologies)[ionto]->getOntologyParser()->getNodesByPosition(ivec)->name << " termPos: " << termPosSucc(ivec) << " termNeg: " << termNegSucc(ivec) << std::endl;
                ++featureNumber;
            }

            if(termPosSucc[ivec] == 0)
                nulla++;
        }
        std::cout << "enrichment for : " << actOnto->getName() << " en count: " << enrichment.count() << " total: " << enrichment.size() << " pocet null:" << nulla << std::endl;

        enrichItems.push_back(enrichment);
        //enrichIDnums.push_back(idNumHash);//refOntologies[ionto]->getOntologyParser()->getNodesByPosition(ibit)->idNum);

    }
    std::vector<boost::dynamic_bitset<> > featureVector(featureNumber);
    for(int iinit = 0; iinit < featureVector.size(); ++iinit)
    {
        //Rcpp::LogicalVector tmp(POSexamples.size() + NEGexamples.size());
        boost::dynamic_bitset<> tmp(POSexamples.size() + NEGexamples.size(), 0);
        featureVector[iinit] = tmp;
    }


    boost::dynamic_bitset<> featureClass(POSexamples.size() + NEGexamples.size());
    std::vector<std::string> featureName(featureNumber);
    this->computeFeatureVectorsbitset(bic, &enrichItems, &featureVector, &featureClass, &featureName);

    for(int ifeat = 0; ifeat < featureVector.size(); ++ifeat)
    {
        enrichPropositionTable.push_back(featureVector[ifeat]);
    }
    //enrichPropositionTable.push_back(featureClass, "Class");
    
	mydataframe newDataTable;
	newDataTable.data = enrichPropositionTable;
	newDataTable.names = featureName;
	newDataTable.predClass = featureClass;
    return newDataTable;
}

Rcpp::List Example::getEnrichTerms()
{
    Rcpp::List enrichmenTerm(refOntologies->size());
    std::vector<std::string> ontologyNames(refOntologies->size());
    for(int ionto = 0; ionto < refOntologies->size(); ++ionto)
    {
        std::vector<double> chistat;
        std::vector<double> pval;
        std::vector<std::string> termID;
        std::vector<std::string> termDesc;
        std::vector<int> nPOS;
        std::vector<int> nNEG;

        int nterm = (*refOntologies)[ionto]->getOntologyParser()->getTermsNumber();
        boost::dynamic_bitset<> mynewclass(this->POSexamples.size() + this->NEGexamples.size(), 0);
        for(int iexample = 0; iexample < this->POSexamples.size(); ++iexample)
        {
            mynewclass[iexample] = 1;
        }

        for(int interm = 0; interm < nterm; ++interm)
        {
            //positive examples
            boost::dynamic_bitset<> mycovered(mynewclass.size(), 0);
            for(int iexample = 0; iexample < this->POSexamples.size(); ++iexample)
            {
                     mycovered[iexample] = (*(this->POSexamples[iexample][ionto]))[interm];//actOnto->getSemanticPatterns();
            }

            //negative examples
            for(int iexample = 0; iexample < this->NEGexamples.size(); ++iexample)
            {
                    mycovered[this->POSexamples.size() + iexample] = (*(this->NEGexamples[iexample][ionto]))[interm];//actOnto->getSemanticPatterns();
            }

            double tp = (mycovered & mynewclass).count();	//(*vector & *classVector).count()
            double fp = (mycovered & ~(mynewclass)).count();	//(*vector & ~(*classVector)).count()
            double posExamples = mynewclass.count();
            double negExamples = (~(mynewclass)).count();

            double res1 = 0;
            if(tp > 0 && posExamples > 0)
                res1 = tp * log2((tp/(tp+fp))/(posExamples/(posExamples+negExamples)));
            double res2 = 0;
            if(fp >0 && negExamples > 0)
                res2 = fp * log2((fp/(tp+fp))/(negExamples/(posExamples+negExamples)));
            chistat.push_back(2*(res1+res2));
            pval.push_back(R::pchisq(2*(res1+res2), 1, false, false));
            nPOS.push_back((int)tp);
            nNEG.push_back((int)fp);
            //Rcpp::Rcout << "pval: " << 2*(res1+res2) << " tp: " << tp << " fp: " << fp << " pos: " << posExamples << " neg: " << negExamples  << std::endl;

            Node *myterm = (*refOntologies)[ionto]->getOntologyParser()->getNodesByPosition(interm);
            termID.push_back(myterm->id);
            termDesc.push_back(myterm->name);
        }

        Rcpp::DataFrame tableEnrich = Rcpp::DataFrame::create(Rcpp::Named("id") =  termID, Rcpp::Named("desc") = termDesc, Rcpp::Named("chistat") = chistat, Rcpp::Named("pval") = pval, Rcpp::Named("npositive") = nPOS, Rcpp::Named("nnegative") = nNEG);
        ontologyNames[ionto] =(*refOntologies)[ionto]->getName();
        enrichmenTerm[ionto] = tableEnrich;

    }
    enrichmenTerm.names() = ontologyNames;
    return enrichmenTerm;
}


mydataframe Example::getEnrichTablev2Sigbitset(int bic, std::vector<Node *> *enrichNodes)
{	
    std::vector<boost::dynamic_bitset<> > enrichPropositionTable;
    std::vector<boost::dynamic_bitset<> > enrichItems;
    //std::vector<boost::unordered_map<unsigned int, bool> > enrichIDnums;

    int featureNumber = 0;
    for(int ionto = 0; ionto < refOntologies->size(); ++ionto)
    {
        int nterm = (*refOntologies)[ionto]->getOntologyParser()->getTermsNumber();
        std::vector<double> pVAL(nterm, 0);

        Ontology *actOnto = (*refOntologies)[ionto];
        //std::vector<boost::dynamic_bitset<> *> patterns = actOnto->getSemanticPatterns();

        //pareto bicluster -> example->ontology->bitset
        this->POSexamples;
        this->NEGexamples;
        
        boost::dynamic_bitset<> mynewclass(this->POSexamples.size() + this->NEGexamples.size(), 0);
        for(int iexample = 0; iexample < this->POSexamples.size(); ++iexample)
        {
			mynewclass[iexample] = 1;
		}

		for(int interm = 0; interm < nterm; ++interm)
		{
			//positive examples
			boost::dynamic_bitset<> mycovered(mynewclass.size(), 0);
			for(int iexample = 0; iexample < this->POSexamples.size(); ++iexample)
			{
				 mycovered[iexample] = (*(this->POSexamples[iexample][ionto]))[interm];//actOnto->getSemanticPatterns();
			}

			//negative examples
			for(int iexample = 0; iexample < this->NEGexamples.size(); ++iexample)
			{
				mycovered[this->POSexamples.size() + iexample] = (*(this->NEGexamples[iexample][ionto]))[interm];//actOnto->getSemanticPatterns();				
			}
			
			double tp = (mycovered & mynewclass).count();	//(*vector & *classVector).count()		
			double fp = (mycovered & ~(mynewclass)).count();	//(*vector & ~(*classVector)).count()
			double posExamples = mynewclass.count();
			double negExamples = (~(mynewclass)).count();

                        //double res1 = tp * log2((tp/(tp+fp))/(posExamples/(posExamples+negExamples)));
                        //double res2 = fp * log2((fp/(tp+fp))/(negExamples/(posExamples+negExamples)));
                        double res1 = 0;
                        if(tp > 0 && posExamples > 0)
                            res1 = tp * log2((tp/(tp+fp))/(posExamples/(posExamples+negExamples)));
                        double res2 = 0;
                        if(fp >0 && negExamples > 0)
                            res2 = fp * log2((fp/(tp+fp))/(negExamples/(posExamples+negExamples)));
			pVAL[interm] = 2*(res1+res2);		
		}
 
        //compute variables
 
        boost::dynamic_bitset<> enrichment(nterm);
//        boost::unordered_map<unsigned int, bool> idNumHash;
        int nulla = 0;
        for(int ivec = 0; ivec < pVAL.size(); ++ivec)
        {
            if(pVAL[ivec] > 6.635)
            {
                enrichment[ivec] = 1;
 //               idNumHash[ivec] = true;
                //refOntologies[ionto].getOntologyParser();
                enrichNodes->push_back((*refOntologies)[ionto]->getOntologyParser()->getNodesByPosition(ivec));
                //Rcpp::Rcout << "ench score: " << pVAL(ivec) << " node: " << (*refOntologies)[ionto]->getOntologyParser()->getNodesByPosition(ivec)->name << " termPos: " << termPosSucc(ivec) << " termNeg: " << termNegSucc(ivec) << std::endl;
                ++featureNumber;
            }

            //if(termPosSucc[ivec] == 0)
            //    nulla++;
        }
        if(this->verbose)
            Rcpp::Rcout << "enrichment for : " << actOnto->getName() << " en count: " << enrichment.count() << " total: " << enrichment.size() << std::endl;

        enrichItems.push_back(enrichment);
        //enrichIDnums.push_back(idNumHash);//refOntologies[ionto]->getOntologyParser()->getNodesByPosition(ibit)->idNum);

    }
    std::vector<boost::dynamic_bitset<> > featureVector(featureNumber);
    for(int iinit = 0; iinit < featureVector.size(); ++iinit)
    {
        //Rcpp::LogicalVector tmp(POSexamples.size() + NEGexamples.size());
        boost::dynamic_bitset<> tmp(POSexamples.size() + NEGexamples.size(), 0);
        featureVector[iinit] = tmp;
    }


    boost::dynamic_bitset<> featureClass(POSexamples.size() + NEGexamples.size());
    std::vector<std::string> featureName(featureNumber);
    this->computeFeatureVectorsbitset(bic, &enrichItems, &featureVector, &featureClass, &featureName);

    for(int ifeat = 0; ifeat < featureVector.size(); ++ifeat)
    {
        enrichPropositionTable.push_back(featureVector[ifeat]);
    }
    //enrichPropositionTable.push_back(featureClass, "Class");
    
	mydataframe newDataTable;
	newDataTable.data = enrichPropositionTable;
	newDataTable.names = featureName;
	newDataTable.predClass = featureClass;
    return newDataTable;
}

Rcpp::DataFrame Example::getEnrichTablev2Ancestors(int bic, std::vector<Node *> *enrichNodes, int minth)
{
    Rcpp::DataFrame enrichPropositionTable;
    std::vector<boost::dynamic_bitset<> > enrichItems;

    int featureNumber = 0;
    for(int ionto = 0; ionto < refOntologies->size(); ++ionto)
    {
        int sizePop = 0;
        int sampleSize = 0;

        int nterm = (*refOntologies)[ionto]->getOntologyParser()->getTermsNumber();
        Rcpp::NumericVector termPosSucc(nterm, 0.0);
        Rcpp::NumericVector termNegSucc(nterm, 0.0);
        Rcpp::NumericVector pVAL(nterm, 0.0);

        Ontology *actOnto = (*refOntologies)[ionto];
        //std::vector<boost::dynamic_bitset<> *> patterns = actOnto->getSemanticPatterns();
        int ppos = 0;
        int nneg = 0;

        //pareto bicluster -> example->ontology->bitset
        this->POSexamples;
        this->NEGexamples;

        //positive examples
        for(int iexample = 0; iexample < this->POSexamples.size(); ++iexample)
        {
            boost::dynamic_bitset<> *patterns = this->POSexamples[iexample][ionto];//actOnto->getSemanticPatterns();


            for(int ibit = 0; ibit < patterns->size(); ++ibit)
            {
                if((*patterns)[ibit] == 1)
                {
                    termPosSucc(ibit) = termPosSucc(ibit) + 1;
                }
            }
            if(patterns->any())
                ++sampleSize;
        }

        //negative examples
        for(int iexample = 0; iexample < this->NEGexamples.size(); ++iexample)
        {
            boost::dynamic_bitset<> *patterns = this->NEGexamples[iexample][ionto];//actOnto->getSemanticPatterns();
            for(int ibit = 0; ibit < patterns->size(); ++ibit)
            {
                if((*patterns)[ibit] == 1)
                {
                    termNegSucc(ibit) = termNegSucc(ibit) + 1;
                }
            }
            if(patterns->any())
                ++sizePop;
        }

        //compute variables
        for(int ivec = 0; ivec < termPosSucc.size(); ++ivec)
        {
            int succSample = termPosSucc(ivec);
            int succPop = termNegSucc(ivec);
            int failurePOP = sizePop-succPop;
            pVAL(ivec) = getHyperGeometricScore(&succSample, &succPop, &failurePOP, &sampleSize);
        }
        //multiple hypothesiss correction - bonferroni
        /* In the case of GO::TermFinder, the value used for the Bonferroni correction is the number of nodes
         *  to which the genes of interest are collectively annotated, excluding those nodes which only have
         *  a single annotation in the background distribution, which a priori cannot be significantly enriched.
        */

        //terms that are collectively annotated
        //arma::uvec indeces = find(termPosSucc == sampleSize);

        //boost::dynamic_bitset<> discovered(nterm);
        int bonferroni = 0;
        for(int ivec = 0; ivec < termPosSucc.size(); ++ivec)
        {
            if(termPosSucc(ivec) >= 1 && termNegSucc(ivec) != 1)
                ++bonferroni;
        }

        std::cout << std::endl;
        double alpha = 0.05;
        std::cout << "alpha: " << alpha << " number of hypothesis: " << bonferroni << " new alpha: " << alpha/(double)bonferroni  << "sample size: " << sampleSize << std::endl;

        boost::dynamic_bitset<> enrichment(nterm);
        int nulla = 0;
        for(int ivec = 0; ivec < pVAL.size(); ++ivec)
        {
            if(pVAL(ivec) <= alpha/(double)bonferroni)
            {
                enrichment[ivec] = 1;
                //refOntologies[ionto].getOntologyParser();
                Node *signode = (*refOntologies)[ionto]->getOntologyParser()->getNodesByPosition(ivec);
                //get all more general nodes
                boost::dynamic_bitset<> moreGeneral = (*refOntologies)[ionto]->getTermBitAncestors(signode);
                enrichment |= moreGeneral;
                //enrichNodes->push_back();
                //std::cout << "ench score: " << pVAL(ivec) << " node: " << (*refOntologies)[ionto]->getOntologyParser()->getNodesByPosition(ivec)->name << " termPos: " << termPosSucc(ivec) << " termNeg: " << termNegSucc(ivec) << std::endl;
                //++featureNumber;
            }

            if(termPosSucc(ivec) == 0)
                nulla++;
        }
        
        boost::dynamic_bitset<> enrichmentAndSize(nterm);
        for(int ibit = 0; ibit < enrichment.size(); ++ibit)
        {
			if(enrichment[ibit] == 1 && (termPosSucc(ibit) + termNegSucc(ibit)) >= minth)	//zde by mela byt kontrola na minimalni pokryti!!! viz termSucc(ivec) >= minth
			{
				enrichNodes->push_back((*refOntologies)[ionto]->getOntologyParser()->getNodesByPosition(ibit));
				enrichmentAndSize[ibit] = 1;
			}
		}
        std::cout << "enrichment for : " << actOnto->getName() << " en count: " << enrichmentAndSize.count() << " total: " << enrichmentAndSize.size() << " pocet null:" << nulla << std::endl;

		featureNumber += enrichmentAndSize.count();
        enrichItems.push_back(enrichmentAndSize);
    }

    /*
     * CO upravit:
     * smazat idNumHash - neni potreba
     * na 817. radku neukladat do enrichnodes, ale pomoci fce getNodesByPosition ziskat node, pak ziskat bitset nedrazenych termu pomoci fce getTermBitAncestors
     * tento vektor |= bit OR operace na enrichment
     * pak znovu iterovat pres cely bitset u ulozit do enrichnodes
     * */


    std::vector<Rcpp::LogicalVector> featureVector(featureNumber);
    for(int iinit = 0; iinit < featureVector.size(); ++iinit)
    {
        Rcpp::LogicalVector tmp(POSexamples.size() + NEGexamples.size());
        featureVector[iinit] = tmp;
    }


    Rcpp::LogicalVector featureClass(POSexamples.size() + NEGexamples.size());
    std::vector<std::string> featureName(featureNumber);
    this->computeFeatureVectors(bic, &enrichItems, &featureVector, &featureClass, &featureName);

    for(int ifeat = 0; ifeat < featureVector.size(); ++ifeat)
    {
        enrichPropositionTable.push_back(featureVector[ifeat], featureName[ifeat]);
    }
    enrichPropositionTable.push_back(featureClass, "Class");

    return enrichPropositionTable;
}



mydataframe Example::getEnrichTablev2Ancestorsbitset(int bic, std::vector<Node *> *enrichNodes, int minth)
{
    std::vector<boost::dynamic_bitset<> > enrichPropositionTable;
    std::vector<boost::dynamic_bitset<> > enrichItems;

    int featureNumber = 0;
    for(int ionto = 0; ionto < refOntologies->size(); ++ionto)
    {
        int sizePop = 0;
        int sampleSize = 0;

        int nterm = (*refOntologies)[ionto]->getOntologyParser()->getTermsNumber();
        std::vector<int> termPosSucc(nterm, 0);
        std::vector<int> termNegSucc(nterm, 0);
        std::vector<double> pVAL(nterm, 0);

        Ontology *actOnto = (*refOntologies)[ionto];
        //std::vector<boost::dynamic_bitset<> *> patterns = actOnto->getSemanticPatterns();
        int ppos = 0;
        int nneg = 0;

        //pareto bicluster -> example->ontology->bitset
        this->POSexamples;
        this->NEGexamples;

        //positive examples
        for(int iexample = 0; iexample < this->POSexamples.size(); ++iexample)
        {
            boost::dynamic_bitset<> *patterns = this->POSexamples[iexample][ionto];//actOnto->getSemanticPatterns();


            for(int ibit = 0; ibit < patterns->size(); ++ibit)
            {
                if((*patterns)[ibit] == 1)
                {
                    termPosSucc[ibit] = termPosSucc[ibit] + 1;
                }
            }
            if(patterns->any())
                ++sampleSize;
        }

        //negative examples
        for(int iexample = 0; iexample < this->NEGexamples.size(); ++iexample)
        {
            boost::dynamic_bitset<> *patterns = this->NEGexamples[iexample][ionto];//actOnto->getSemanticPatterns();
            for(int ibit = 0; ibit < patterns->size(); ++ibit)
            {
                if((*patterns)[ibit] == 1)
                {
                    termNegSucc[ibit] = termNegSucc[ibit] + 1;
                }
            }
            if(patterns->any())
                ++sizePop;
        }

        //compute variables
        for(int ivec = 0; ivec < termPosSucc.size(); ++ivec)
        {
            int succSample = termPosSucc[ivec];
            int succPop = termNegSucc[ivec];
            int failurePOP = sizePop-succPop;
            pVAL[ivec] = getHyperGeometricScore(&succSample, &succPop, &failurePOP, &sampleSize);
        }
        //multiple hypothesiss correction - bonferroni
        /* In the case of GO::TermFinder, the value used for the Bonferroni correction is the number of nodes
         *  to which the genes of interest are collectively annotated, excluding those nodes which only have
         *  a single annotation in the background distribution, which a priori cannot be significantly enriched.
        */

        //terms that are collectively annotated
        //arma::uvec indeces = find(termPosSucc == sampleSize);

        //boost::dynamic_bitset<> discovered(nterm);
        int bonferroni = 0;
        for(int ivec = 0; ivec < termPosSucc.size(); ++ivec)
        {
            if(termPosSucc[ivec] >= 1 && termNegSucc[ivec] != 1)
                ++bonferroni;
        }

        std::cout << std::endl;
        double alpha = 0.05;
        std::cout << "alpha: " << alpha << " number of hypothesis: " << bonferroni << " new alpha: " << alpha/(double)bonferroni  << "sample size: " << sampleSize << std::endl;

        boost::dynamic_bitset<> enrichment(nterm);
        int nulla = 0;
        for(int ivec = 0; ivec < pVAL.size(); ++ivec)
        {
            if(pVAL[ivec] <= alpha/(double)bonferroni)
            {
                enrichment[ivec] = 1;
                //refOntologies[ionto].getOntologyParser();
                Node *signode = (*refOntologies)[ionto]->getOntologyParser()->getNodesByPosition(ivec);
                //get all more general nodes
                boost::dynamic_bitset<> moreGeneral = (*refOntologies)[ionto]->getTermBitAncestors(signode);
                enrichment |= moreGeneral;
                //enrichNodes->push_back();
                //std::cout << "ench score: " << pVAL(ivec) << " node: " << (*refOntologies)[ionto]->getOntologyParser()->getNodesByPosition(ivec)->name << " termPos: " << termPosSucc(ivec) << " termNeg: " << termNegSucc(ivec) << std::endl;
                //++featureNumber;
            }

            if(termPosSucc[ivec] == 0)
                nulla++;
        }
        
        boost::dynamic_bitset<> enrichmentAndSize(nterm);
        for(int ibit = 0; ibit < enrichment.size(); ++ibit)
        {
			if(enrichment[ibit] == 1 && (termPosSucc[ibit] + termNegSucc[ibit]) >= minth)	//zde by mela byt kontrola na minimalni pokryti!!! viz termSucc(ivec) >= minth
			{
				enrichNodes->push_back((*refOntologies)[ionto]->getOntologyParser()->getNodesByPosition(ibit));
				enrichmentAndSize[ibit] = 1;
			}
		}
        std::cout << "enrichment for : " << actOnto->getName() << " en count: " << enrichmentAndSize.count() << " total: " << enrichmentAndSize.size() << " pocet null:" << nulla << std::endl;

		featureNumber += enrichmentAndSize.count();
        enrichItems.push_back(enrichmentAndSize);
    }

    /*
     * CO upravit:
     * smazat idNumHash - neni potreba
     * na 817. radku neukladat do enrichnodes, ale pomoci fce getNodesByPosition ziskat node, pak ziskat bitset nedrazenych termu pomoci fce getTermBitAncestors
     * tento vektor |= bit OR operace na enrichment
     * pak znovu iterovat pres cely bitset u ulozit do enrichnodes
     * */


    std::vector<boost::dynamic_bitset<> > featureVector(featureNumber);
    for(int iinit = 0; iinit < featureVector.size(); ++iinit)
    {
        //Rcpp::LogicalVector tmp(POSexamples.size() + NEGexamples.size());
        boost::dynamic_bitset<> tmp(POSexamples.size() + NEGexamples.size(), 0);
        featureVector[iinit] = tmp;
    }


    boost::dynamic_bitset<> featureClass(POSexamples.size() + NEGexamples.size());
    std::vector<std::string> featureName(featureNumber);
    this->computeFeatureVectorsbitset(bic, &enrichItems, &featureVector, &featureClass, &featureName);

    for(int ifeat = 0; ifeat < featureVector.size(); ++ifeat)
    {
        enrichPropositionTable.push_back(featureVector[ifeat]);
    }
    //enrichPropositionTable.push_back(featureClass, "Class");
    
	mydataframe newDataTable;
	newDataTable.data = enrichPropositionTable;
	newDataTable.names = featureName;
	newDataTable.predClass = featureClass;
    
    return newDataTable;
}


mydataframe Example::getEnrichTablev2AncestorsSigbitset(int bic, std::vector<Node *> *enrichNodes, int minth)
{
    std::vector<boost::dynamic_bitset<> > enrichPropositionTable;
    std::vector<boost::dynamic_bitset<> > enrichItems;

    int featureNumber = 0;
    for(int ionto = 0; ionto < refOntologies->size(); ++ionto)
    {
        int sizePop = 0;
        int sampleSize = 0;

        int nterm = (*refOntologies)[ionto]->getOntologyParser()->getTermsNumber();
        std::vector<int> termPosSucc(nterm, 0);
        std::vector<int> termNegSucc(nterm, 0);
        std::vector<double> pVAL(nterm, 0);

        Ontology *actOnto = (*refOntologies)[ionto];
        //std::vector<boost::dynamic_bitset<> *> patterns = actOnto->getSemanticPatterns();

        //pareto bicluster -> example->ontology->bitset
        this->POSexamples;
        this->NEGexamples;
        
		boost::dynamic_bitset<> mynewclass(this->POSexamples.size() + this->NEGexamples.size(), 0);
        for(int iexample = 0; iexample < this->POSexamples.size(); ++iexample)
        {
			mynewclass[iexample] = 1;
		}

		for(int interm = 0; interm < nterm; ++interm)
		{
			//positive examples
			boost::dynamic_bitset<> mycovered(mynewclass.size(), 0);
			for(int iexample = 0; iexample < this->POSexamples.size(); ++iexample)
			{
				 mycovered[iexample] = (*(this->POSexamples[iexample][ionto]))[interm];//actOnto->getSemanticPatterns();
			}

			//negative examples
			for(int iexample = 0; iexample < this->NEGexamples.size(); ++iexample)
			{
				mycovered[this->POSexamples.size() + iexample] = (*(this->NEGexamples[iexample][ionto]))[interm];//actOnto->getSemanticPatterns();				
			}
			
			double tp = (mycovered & mynewclass).count();	//(*vector & *classVector).count()		
			double fp = (mycovered & ~(mynewclass)).count();	//(*vector & ~(*classVector)).count()
			double posExamples = mynewclass.count();
			double negExamples = (~(mynewclass)).count();

                        //double res1 = tp * log2((tp/(tp+fp))/(posExamples/(posExamples+negExamples)));
                        //double res2 = fp * log2((fp/(tp+fp))/(negExamples/(posExamples+negExamples)));
                        double res1 = 0;
                        if(tp > 0 && posExamples > 0)
                            res1 = tp * log2((tp/(tp+fp))/(posExamples/(posExamples+negExamples)));
                        double res2 = 0;
                        if(fp >0 && negExamples > 0)
                            res2 = fp * log2((fp/(tp+fp))/(negExamples/(posExamples+negExamples)));
			pVAL[interm] = 2*(res1+res2);		
		}

        //positive examples
        for(int iexample = 0; iexample < this->POSexamples.size(); ++iexample)
        {
            boost::dynamic_bitset<> *patterns = this->POSexamples[iexample][ionto];//actOnto->getSemanticPatterns();


            for(int ibit = 0; ibit < patterns->size(); ++ibit)
            {
                if((*patterns)[ibit] == 1)
                {
                    termPosSucc[ibit] = termPosSucc[ibit] + 1;
                }
            }
            if(patterns->any())
                ++sampleSize;
        }

        //negative examples
        for(int iexample = 0; iexample < this->NEGexamples.size(); ++iexample)
        {
            boost::dynamic_bitset<> *patterns = this->NEGexamples[iexample][ionto];//actOnto->getSemanticPatterns();
            for(int ibit = 0; ibit < patterns->size(); ++ibit)
            {
                if((*patterns)[ibit] == 1)
                {
                    termNegSucc[ibit] = termNegSucc[ibit] + 1;
                }
            }
            if(patterns->any())
                ++sizePop;
        }

        Rcpp::Rcout << std::endl;
        double alpha = 0.05;
        //std::cout << "alpha: " << alpha << " number of hypothesis: " << bonferroni << " new alpha: " << alpha/(double)bonferroni  << "sample size: " << sampleSize << std::endl;

        boost::dynamic_bitset<> enrichment(nterm);
        int nulla = 0;
        for(int ivec = 0; ivec < pVAL.size(); ++ivec)
        {
            if(pVAL[ivec] > 6.635)
            {
                enrichment[ivec] = 1;
                //refOntologies[ionto].getOntologyParser();
                Node *signode = (*refOntologies)[ionto]->getOntologyParser()->getNodesByPosition(ivec);
                //get all more general nodes
                boost::dynamic_bitset<> moreGeneral = (*refOntologies)[ionto]->getTermBitAncestors(signode);
                enrichment |= moreGeneral;
                //enrichNodes->push_back();
                //std::cout << "ench score: " << pVAL(ivec) << " node: " << (*refOntologies)[ionto]->getOntologyParser()->getNodesByPosition(ivec)->name << " termPos: " << termPosSucc(ivec) << " termNeg: " << termNegSucc(ivec) << std::endl;
                //++featureNumber;
            }

            if(termPosSucc[ivec] == 0)
                nulla++;
        }
        
        boost::dynamic_bitset<> enrichmentAndSize(nterm);
        for(int ibit = 0; ibit < enrichment.size(); ++ibit)
        {
			if(enrichment[ibit] == 1 && (termPosSucc[ibit] + termNegSucc[ibit]) >= minth)	//zde by mela byt kontrola na minimalni pokryti!!! viz termSucc(ivec) >= minth
			{
				enrichNodes->push_back((*refOntologies)[ionto]->getOntologyParser()->getNodesByPosition(ibit));
				enrichmentAndSize[ibit] = 1;
			}
		}
        if(this->verbose)
            Rcpp::Rcout << "enrichment for : " << actOnto->getName() << " en count: " << enrichmentAndSize.count() << " total: " << enrichmentAndSize.size() << " pocet null:" << nulla << std::endl;

		featureNumber += enrichmentAndSize.count();
        enrichItems.push_back(enrichmentAndSize);
    }

    /*
     * CO upravit:
     * smazat idNumHash - neni potreba
     * na 817. radku neukladat do enrichnodes, ale pomoci fce getNodesByPosition ziskat node, pak ziskat bitset nedrazenych termu pomoci fce getTermBitAncestors
     * tento vektor |= bit OR operace na enrichment
     * pak znovu iterovat pres cely bitset u ulozit do enrichnodes
     * */


    std::vector<boost::dynamic_bitset<> > featureVector(featureNumber);
    for(int iinit = 0; iinit < featureVector.size(); ++iinit)
    {
        //Rcpp::LogicalVector tmp(POSexamples.size() + NEGexamples.size());
        boost::dynamic_bitset<> tmp(POSexamples.size() + NEGexamples.size(), 0);
        featureVector[iinit] = tmp;
    }


    boost::dynamic_bitset<> featureClass(POSexamples.size() + NEGexamples.size());
    std::vector<std::string> featureName(featureNumber);
    this->computeFeatureVectorsbitset(bic, &enrichItems, &featureVector, &featureClass, &featureName);

    for(int ifeat = 0; ifeat < featureVector.size(); ++ifeat)
    {
        enrichPropositionTable.push_back(featureVector[ifeat]);
    }
    //enrichPropositionTable.push_back(featureClass, "Class");
    
	mydataframe newDataTable;
	newDataTable.data = enrichPropositionTable;
	newDataTable.names = featureName;
	newDataTable.predClass = featureClass;
    
    return newDataTable;
}


Rcpp::DataFrame Example::getEnrichTablev3TH(int bic, std::vector<Node *> *enrichNodes, int minth)
{
    Rcpp::DataFrame enrichPropositionTable;
    std::vector<boost::dynamic_bitset<> > enrichItems;
    std::vector<boost::unordered_map<unsigned int, bool> > enrichIDnums;

    int featureNumber = 0;
    for(int ionto = 0; ionto < refOntologies->size(); ++ionto)
    {
        int sizePop = 0;
        int sampleSize = 0;

        int nterm = (*refOntologies)[ionto]->getOntologyParser()->getTermsNumber();
        Rcpp::NumericVector termSucc(nterm, 0.0);

        Ontology *actOnto = (*refOntologies)[ionto];
        //std::vector<boost::dynamic_bitset<> *> patterns = actOnto->getSemanticPatterns();

        //pareto bicluster -> example->ontology->bitset
        this->POSexamples;
        this->NEGexamples;

        //positive examples
        for(int iexample = 0; iexample < this->POSexamples.size(); ++iexample)
        {
            boost::dynamic_bitset<> *patterns = this->POSexamples[iexample][ionto];//actOnto->getSemanticPatterns();
            for(int ibit = 0; ibit < patterns->size(); ++ibit)
            {
                if((*patterns)[ibit] == 1)
                {
                    termSucc(ibit) = termSucc(ibit) + 1;
                }
            }
        }

        //negative examples
        for(int iexample = 0; iexample < this->NEGexamples.size(); ++iexample)
        {
            boost::dynamic_bitset<> *patterns = this->NEGexamples[iexample][ionto];//actOnto->getSemanticPatterns();
            for(int ibit = 0; ibit < patterns->size(); ++ibit)
            {
                if((*patterns)[ibit] == 1)
                {
                    termSucc(ibit) = termSucc(ibit) + 1;
                }
            }
        }

        boost::dynamic_bitset<> enrichment(nterm);
        boost::unordered_map<unsigned int, bool> idNumHash;
        int nulla = 0;
        for(int ivec = 0; ivec < termSucc.size(); ++ivec)
        {
            if(termSucc(ivec) >= minth)
            {
                enrichment[ivec] = 1;
                idNumHash[ivec] = true;
                //refOntologies[ionto].getOntologyParser();
                enrichNodes->push_back((*refOntologies)[ionto]->getOntologyParser()->getNodesByPosition(ivec));
                //std::cout << "ench score: " << pVAL(ivec) << " node: " << (*refOntologies)[ionto]->getOntologyParser()->getNodesByPosition(ivec)->name << " termPos: " << termPosSucc(ivec) << " termNeg: " << termNegSucc(ivec) << std::endl;
                ++featureNumber;
            }

            if(termSucc(ivec) == 0)
                nulla++;
        }
        std::cout << "enrichment for : " << actOnto->getName() << " en count: " << enrichment.count() << " total: " << enrichment.size() << " pocet null:" << nulla << std::endl;

        enrichItems.push_back(enrichment);
        enrichIDnums.push_back(idNumHash);//refOntologies[ionto]->getOntologyParser()->getNodesByPosition(ibit)->idNum);

    }
    std::vector<Rcpp::LogicalVector> featureVector(featureNumber);
    for(int iinit = 0; iinit < featureVector.size(); ++iinit)
    {
        Rcpp::LogicalVector tmp(POSexamples.size() + NEGexamples.size());
        featureVector[iinit] = tmp;
    }


    Rcpp::LogicalVector featureClass(POSexamples.size() + NEGexamples.size());
    std::vector<std::string> featureName(featureNumber);
    this->computeFeatureVectors(bic, &enrichItems, &featureVector, &featureClass, &featureName);

    for(int ifeat = 0; ifeat < featureVector.size(); ++ifeat)
    {
        enrichPropositionTable.push_back(featureVector[ifeat], featureName[ifeat]);
    }
    enrichPropositionTable.push_back(featureClass, "Class");

    return enrichPropositionTable;
}

mydataframe Example::getEnrichTablev3THbitset(int bic, std::vector<Node *> *enrichNodes, int minth)
{
    std::vector<boost::dynamic_bitset<> > enrichPropositionTable;
    std::vector<boost::dynamic_bitset<> > enrichItems;
    //std::vector<boost::unordered_map<unsigned int, bool> > enrichIDnums;

    int featureNumber = 0;
    for(int ionto = 0; ionto < refOntologies->size(); ++ionto)
    {
        int nterm = (*refOntologies)[ionto]->getOntologyParser()->getTermsNumber();
        std::vector<int> termSucc(nterm, 0.0);

        Ontology *actOnto = (*refOntologies)[ionto];
        //std::vector<boost::dynamic_bitset<> *> patterns = actOnto->getSemanticPatterns();

        //pareto bicluster -> example->ontology->bitset
        this->POSexamples;
        this->NEGexamples;

        //positive examples
        for(int iexample = 0; iexample < this->POSexamples.size(); ++iexample)
        {
            boost::dynamic_bitset<> *patterns = this->POSexamples[iexample][ionto];//actOnto->getSemanticPatterns();
            for(int ibit = 0; ibit < patterns->size(); ++ibit)
            {
                if((*patterns)[ibit] == 1)
                {
                    termSucc[ibit] = termSucc[ibit] + 1;
                }
            }
        }

        //negative examples
        for(int iexample = 0; iexample < this->NEGexamples.size(); ++iexample)
        {
            boost::dynamic_bitset<> *patterns = this->NEGexamples[iexample][ionto];//actOnto->getSemanticPatterns();
            for(int ibit = 0; ibit < patterns->size(); ++ibit)
            {
                if((*patterns)[ibit] == 1)
                {
                    termSucc[ibit] = termSucc[ibit] + 1;
                }
            }
        }

        boost::dynamic_bitset<> enrichment(nterm);
        //boost::unordered_map<unsigned int, bool> idNumHash;
        int nulla = 0;
        for(int ivec = 0; ivec < termSucc.size(); ++ivec)
        {
            if(termSucc[ivec] >= minth)
            {
                enrichment[ivec] = 1;
                //idNumHash[ivec] = true;
                //refOntologies[ionto].getOntologyParser();
                enrichNodes->push_back((*refOntologies)[ionto]->getOntologyParser()->getNodesByPosition(ivec));
                //std::cout << "ench score: " << pVAL(ivec) << " node: " << (*refOntologies)[ionto]->getOntologyParser()->getNodesByPosition(ivec)->name << " termPos: " << termPosSucc(ivec) << " termNeg: " << termNegSucc(ivec) << std::endl;
                ++featureNumber;
            }

            if(termSucc[ivec] == 0)
                nulla++;
        }
        if(this->verbose)
            Rcpp::Rcout << "enrichment for : " << actOnto->getName() << " en count: " << enrichment.count() << " total: " << enrichment.size() << " pocet null:" << nulla << std::endl;

        enrichItems.push_back(enrichment);
        //enrichIDnums.push_back(idNumHash);//refOntologies[ionto]->getOntologyParser()->getNodesByPosition(ibit)->idNum);

    }
    
    std::vector<boost::dynamic_bitset<> > featureVector(featureNumber);
    for(int iinit = 0; iinit < featureVector.size(); ++iinit)
    {
        //Rcpp::LogicalVector tmp(POSexamples.size() + NEGexamples.size());
        boost::dynamic_bitset<> tmp(POSexamples.size() + NEGexamples.size(), 0);
        featureVector[iinit] = tmp;
    }


    boost::dynamic_bitset<> featureClass(POSexamples.size() + NEGexamples.size());
    std::vector<std::string> featureName(featureNumber);
    this->computeFeatureVectorsbitset(bic, &enrichItems, &featureVector, &featureClass, &featureName);

    for(int ifeat = 0; ifeat < featureVector.size(); ++ifeat)
    {
        enrichPropositionTable.push_back(featureVector[ifeat]);
    }
    //enrichPropositionTable.push_back(featureClass, "Class");
    
	mydataframe newDataTable;
	newDataTable.data = enrichPropositionTable;
	newDataTable.names = featureName;
	newDataTable.predClass = featureClass;
	
    return newDataTable;
}


void Example::getEnrichTablev3THNegative(int bic, std::vector<Node *> *enrichNodes, int minth, Rcpp::DataFrame *enrichPropositionTable)
{
    //Rcpp::DataFrame enrichPropositionTable;
    //erase class
    enrichPropositionTable->erase(enrichPropositionTable->length()-1);
    std::vector<boost::dynamic_bitset<> > enrichItems;
    std::vector<boost::unordered_map<unsigned int, bool> > enrichIDnums;

    int featureNumber = 0;
    for(int ionto = 0; ionto < refOntologies->size(); ++ionto)
    {
        int sizePop = 0;
        int sampleSize = 0;

        int nterm = (*refOntologies)[ionto]->getOntologyParser()->getTermsNumber();
        Rcpp::NumericVector termNegSucc(nterm, 0.0);

        Ontology *actOnto = (*refOntologies)[ionto];
        //std::vector<boost::dynamic_bitset<> *> patterns = actOnto->getSemanticPatterns();

        //pareto bicluster -> example->ontology->bitset
        this->POSexamples;
        this->NEGexamples;

        //positive examples
        for(int iexample = 0; iexample < this->POSexamples.size(); ++iexample)
        {
            boost::dynamic_bitset<> *patterns = this->POSexamples[iexample][ionto];//actOnto->getSemanticPatterns();
            for(int ibit = 0; ibit < patterns->size(); ++ibit)
            {
                //if((*patterns)[ibit] == 0)
                if((*patterns)[ibit] == 1)
                {
                    termNegSucc(ibit) = termNegSucc(ibit) + 1;
                }
            }
        }

        //negative examples
        for(int iexample = 0; iexample < this->NEGexamples.size(); ++iexample)
        {
            boost::dynamic_bitset<> *patterns = this->NEGexamples[iexample][ionto];//actOnto->getSemanticPatterns();
            for(int ibit = 0; ibit < patterns->size(); ++ibit)
            {
                /*if((*patterns)[ibit] == 0)
                {
                    termNegSucc(ibit) = termNegSucc(ibit) - 1;
                }*/
            }
        }

        boost::dynamic_bitset<> enrichment(nterm);
        boost::unordered_map<unsigned int, bool> idNumHash;
        int nulla = 0;
        for(int ivec = 0; ivec < termNegSucc.size(); ++ivec)
        {
            if(termNegSucc(ivec) >= minth)
            {
                enrichment[ivec] = 1;
                idNumHash[ivec] = true;
                //refOntologies[ionto].getOntologyParser();
                enrichNodes->push_back((*refOntologies)[ionto]->getOntologyParser()->getNodesByPosition_bottomUP(ivec));
                //std::cout << "ench score: " << pVAL(ivec) << " node: " << (*refOntologies)[ionto]->getOntologyParser()->getNodesByPosition(ivec)->name << " termPos: " << termPosSucc(ivec) << " termNeg: " << termNegSucc(ivec) << std::endl;
                ++featureNumber;
            }

            if(termNegSucc(ivec) == 0)
                nulla++;
        }
        std::cout << "negative enrichment for : " << actOnto->getName() << " en count: " << enrichment.count() << " total: " << enrichment.size() << " pocet null:" << nulla << std::endl;

        enrichItems.push_back(enrichment);
        enrichIDnums.push_back(idNumHash);//refOntologies[ionto]->getOntologyParser()->getNodesByPosition(ibit)->idNum);

    }
    std::vector<Rcpp::LogicalVector> featureVector(featureNumber);
    for(int iinit = 0; iinit < featureVector.size(); ++iinit)
    {
        Rcpp::LogicalVector tmp(POSexamples.size() + NEGexamples.size());
        featureVector[iinit] = tmp;
    }


    Rcpp::LogicalVector featureClass(POSexamples.size() + NEGexamples.size());
    std::vector<std::string> featureName(featureNumber);
    this->computeFeatureVectors(bic, &enrichItems, &featureVector, &featureClass, &featureName);

    for(int ifeat = 0; ifeat < featureVector.size(); ++ifeat)
    {
        enrichPropositionTable->push_back(featureVector[ifeat], featureName[ifeat]);
    }
    enrichPropositionTable->push_back(featureClass, "Class");
}

boost::dynamic_bitset<> Example::getEnrichClass(Rcpp::DataFrame *enrichPropositionTable)
{
    Rcpp::NumericVector myclass = (*enrichPropositionTable)["Class"];
    boost::dynamic_bitset<> newclass(myclass.size(), 0);
    for(int i = 0; i < newclass.size(); ++i)
    {
        if(myclass[i] == 1)
            newclass[i] = 1;
    }
    return newclass;
}

void Example::computeFeatureVectors(int bic, std::vector<boost::dynamic_bitset<> > *enrichItems, std::vector<Rcpp::LogicalVector> *featureVector, Rcpp::LogicalVector *featureClass, std::vector<std::string> *featureName)
{
    //->feature -> example
    //std::vector<Rcpp::LogicalVector> featureVector(featureNumber);
    //featureVector.resize(featureNumber);
    //Rcpp::LogicalVector featureClass(POSexamples[bic].size() + NEGexamples.size());

    for(int ionto = 0; ionto < refOntologies->size(); ++ionto)
    {
        for(int iench = 0; iench < (*enrichItems)[ionto].size(); ++iench)
        {
            //enriched item
            if((*enrichItems)[ionto][iench] == 1)
            {
                featureName->push_back((*refOntologies)[ionto]->getOntologyParser()->getNodesByPosition(iench)->name);
            }
        }
    }
    //bicluster 0 - positive examples
    for(int iexpos = 0; iexpos < POSexamples.size(); ++iexpos)
    {
        int ivect = 0;
        for(int ionto = 0; ionto < POSexamples[iexpos].size(); ++ionto)
        {
            for(int iench = 0; iench < (*enrichItems)[ionto].size(); ++iench)
            {
                //enriched item
                if((*enrichItems)[ionto][iench] == 1)
                {
                    //return;
                    //std::cout << "iech: " << iench << " ivect: " << ivect << " possize: " << featureVector.size() <<;

                    if(POSexamples[iexpos][ionto]->test(iench))
                        (*featureVector)[ivect][iexpos] = true;
                    else
                        (*featureVector)[ivect][iexpos] = false;

                   ++ivect;
                }
            }
        }
        (*featureClass)[iexpos] = true;
    }    
    //bicluster 0 - negative examples
    for(int iexpos = 0; iexpos < NEGexamples.size(); ++iexpos)
    {
        int ivect = 0;
        for(int ionto = 0; ionto < NEGexamples[iexpos].size(); ++ionto)
        {
            for(int iench = 0; iench < (*enrichItems)[ionto].size(); ++iench)
            {
                //enriched item
                if((*enrichItems)[ionto][iench] == 1)
                {
                    if(NEGexamples[iexpos][ionto]->test(iench))
                        (*featureVector)[ivect][iexpos+POSexamples.size()] = true;
                    else
                        (*featureVector)[ivect][iexpos+POSexamples.size()] = false;
                    ++ivect;
                }                
            }            
        }
        (*featureClass)[iexpos + POSexamples.size()] = false;
    }
}

void Example::computeFeatureVectorsbitset(int bic, std::vector<boost::dynamic_bitset<> > *enrichItems, std::vector<boost::dynamic_bitset<> > *featureVector, boost::dynamic_bitset<> *featureClass, std::vector<std::string> *featureName)
{
    //->feature -> example
    //std::vector<Rcpp::LogicalVector> featureVector(featureNumber);
    //featureVector.resize(featureNumber);
    //Rcpp::LogicalVector featureClass(POSexamples[bic].size() + NEGexamples.size());

    for(int ionto = 0; ionto < refOntologies->size(); ++ionto)
    {
        for(int iench = 0; iench < (*enrichItems)[ionto].size(); ++iench)
        {
            //enriched item
            if((*enrichItems)[ionto][iench] == 1)
            {
                featureName->push_back((*refOntologies)[ionto]->getOntologyParser()->getNodesByPosition(iench)->name);
            }
        }
    }
    //bicluster 0 - positive examples
    for(int iexpos = 0; iexpos < POSexamples.size(); ++iexpos)
    {
        int ivect = 0;
        for(int ionto = 0; ionto < POSexamples[iexpos].size(); ++ionto)
        {
            for(int iench = 0; iench < (*enrichItems)[ionto].size(); ++iench)
            {
                //enriched item
                if((*enrichItems)[ionto][iench] == 1)
                {
                    //return;
                    //std::cout << "iech: " << iench << " ivect: " << ivect << " possize: " << featureVector.size() <<;

                    if(POSexamples[iexpos][ionto]->test(iench))
                        (*featureVector)[ivect][iexpos] = true;
                    else
                        (*featureVector)[ivect][iexpos] = false;

                   ++ivect;
                }
            }
        }
        (*featureClass)[iexpos] = true;
    }    
    //bicluster 0 - negative examples
    for(int iexpos = 0; iexpos < NEGexamples.size(); ++iexpos)
    {
        int ivect = 0;
        for(int ionto = 0; ionto < NEGexamples[iexpos].size(); ++ionto)
        {
            for(int iench = 0; iench < (*enrichItems)[ionto].size(); ++iench)
            {
                //enriched item
                if((*enrichItems)[ionto][iench] == 1)
                {
                    if(NEGexamples[iexpos][ionto]->test(iench))
                        (*featureVector)[ivect][iexpos+POSexamples.size()] = true;
                    else
                        (*featureVector)[ivect][iexpos+POSexamples.size()] = false;
                    ++ivect;
                }                
            }            
        }
        (*featureClass)[iexpos + POSexamples.size()] = false;
    }
}

double Example::getHyperGeometricScore(int *successSample, int *successPop, int *failurePop, int *sampleSize)
{
    //http://users.unimi.it/marray/2007/material/day4/Lecture7.pdf
    int log = 0;
    int lower = 0;
    return R::phyper(*successSample - 1, *successPop, *failurePop, *sampleSize, lower, log );
}

std::vector<bottomFeature> Example::initBottomFeatures(std::vector<Node *> *enrichNodes, Rcpp::DataFrame *enrichPropositionTable, int onlyPosFeatures)
{
    std::vector<bottomFeature> bottomRules(enrichNodes->size());
    int bitsetSize = enrichNodes->size();
    if(!onlyPosFeatures)
            bitsetSize = 2*bitsetSize;

    boost::unordered_map<Node*, int> nodePosition;
    for(int i = 0; i < bottomRules.size(); ++i)
    {
        bottomFeature tmp;        
        Rcpp::NumericVector vctCol = (*enrichPropositionTable)[i];
        //set cover bitset
        boost::dynamic_bitset<> coverbitset(vctCol.size(), 0);
        for(int ibit = 0; ibit < vctCol.size(); ++ibit)
        {
            if(vctCol[ibit] == 1)
                coverbitset[ibit] = 1;
        }
        tmp.exampleCovered = coverbitset;
        tmp.IDind = i;
        tmp.posFeature = 1;
        tmp.actScore = 0;
        tmp.pathBestScore = 0;
        tmp.noderef = (*enrichNodes)[i];
        tmp.allChilds = boost::dynamic_bitset<>(bitsetSize, 0);
        tmp.allParents = boost::dynamic_bitset<>(bitsetSize, 0);
        tmp.allSpecific = boost::dynamic_bitset<>(bitsetSize, 0);
        tmp.allGeneral = boost::dynamic_bitset<>(bitsetSize, 0);
        tmp.level = 0;

        nodePosition[tmp.noderef] = i;

        bottomRules[i] = tmp;
    }

    boost::unordered_map<Node*, int> nodePositionNEG;
    if(!onlyPosFeatures)
    {
        int origVectSize = bottomRules.size();
        bottomRules.resize(2*origVectSize);
        for(int i = 0; i < origVectSize; ++i)
        {
            bottomFeature tmp;
            Rcpp::NumericVector vctCol = (*enrichPropositionTable)[i];
            //set cover bitset
            boost::dynamic_bitset<> coverbitset(vctCol.size(), 0);
            for(int ibit = 0; ibit < vctCol.size(); ++ibit)
            {
                //negative
                if(vctCol[ibit] == 0)
                    coverbitset[ibit] = 1;
            }
            tmp.exampleCovered = coverbitset;
            tmp.IDind = i+origVectSize;
            tmp.posFeature = 0;
            tmp.actScore = 0;
            tmp.pathBestScore = 0;
            tmp.noderef = (*enrichNodes)[i];
            tmp.allChilds = boost::dynamic_bitset<>(bitsetSize, 0);
            tmp.allParents = boost::dynamic_bitset<>(bitsetSize, 0);
            tmp.allSpecific = boost::dynamic_bitset<>(bitsetSize, 0);
            tmp.allGeneral = boost::dynamic_bitset<>(bitsetSize, 0);
            tmp.level = 0;

            nodePositionNEG[tmp.noderef] = i+origVectSize;

            bottomRules[i+origVectSize] = tmp;
        }
    }

    //reduce unused nodes in heirarchies
    for(int iench = 0; iench < enrichNodes->size(); ++iench)
    {
        std::list<Node*> OPEN;
        OPEN.push_back((*enrichNodes)[iench]);
        Node *basedNode = OPEN.front();
        boost::unordered_map<Node*, bool> CLOSED;
        while(!OPEN.empty())
        {
            Node *actNode = OPEN.front();
            OPEN.pop_front();
            std::vector<edge> childs = actNode->relationship;
            for(int ichild = 0; ichild < childs.size(); ++ichild)
            {
                //is not a bottom and is not in CLOSED
                if(childs[ichild].end != NULL && CLOSED.find(childs[ichild].end) == CLOSED.end())
                {
                    //is the enriched node
                    if(nodePosition.find(childs[ichild].end) != nodePosition.end())
                    {
                        //relation between child and parent
                        bottomRules[nodePosition[childs[ichild].end]].allParents[nodePosition[basedNode]] = 1;
                        //relation between parent and child
                        bottomRules[nodePosition[basedNode]].allChilds[nodePosition[childs[ichild].end]] = 1;
                    }
                    else    //find another node
                    {
                        OPEN.push_back(childs[ichild].end);
                    }
                }
            }
            CLOSED[actNode] = true;
        }
    }

    if(!onlyPosFeatures)
    {
        //reduce unused nodes in heirarchies
        for(int iench = 0; iench < enrichNodes->size(); ++iench)
        {
            std::list<Node*> OPEN;
            OPEN.push_back((*enrichNodes)[iench]);
            Node *basedNode = OPEN.front();
            boost::unordered_map<Node*, bool> CLOSED;
            while(!OPEN.empty())
            {
                Node *actNode = OPEN.front();
                OPEN.pop_front();
                std::vector<edge> childs = actNode->relationship;
                for(int ichild = 0; ichild < childs.size(); ++ichild)
                {
                    //is not a bottom and is not in CLOSED
                    if(childs[ichild].end != NULL && CLOSED.find(childs[ichild].end) == CLOSED.end())
                    {
                        //is the enriched node
                        if(nodePositionNEG.find(childs[ichild].end) != nodePositionNEG.end())
                        {
                            //relation between child and parent
                            bottomRules[nodePositionNEG[childs[ichild].end]].allParents[nodePositionNEG[basedNode]] = 1;
                            //relation between parent and child
                            bottomRules[nodePositionNEG[basedNode]].allChilds[nodePositionNEG[childs[ichild].end]] = 1;
                        }
                        else    //find another node
                        {
                            OPEN.push_back(childs[ichild].end);
                        }
                    }
                }
                CLOSED[actNode] = true;
            }
        }
    }

    //complete all specific and general nodes
    for(int ibottom = 0; ibottom < bottomRules.size(); ++ibottom)
    {
        bottomRules[ibottom].allSpecific = completAllChilds(&(bottomRules[ibottom]), &bottomRules);
        bottomRules[ibottom].allGeneral = completAllParents(&(bottomRules[ibottom]), &bottomRules);
    }
    
        //find all roots
    std::list<int> OPEN;
    for(int ibottom = 0; ibottom < bottomRules.size(); ++ibottom)
    {
		if(bottomRules[ibottom].allGeneral.none())
		{
			OPEN.push_back(ibottom);			
		}
	}
	
	boost::dynamic_bitset<> CLOSED(bitsetSize, 0);
	while(!OPEN.empty())
	{
		int ibottom = OPEN.front();
		OPEN.pop_front();
		size_t index = bottomRules[ibottom].allSpecific.find_first();
		
		while(index != boost::dynamic_bitset<>::npos)
		{
			if(CLOSED[index] == 0)
			{
				
				OPEN.push_back(index);
                                bottomRules[index].level = bottomRules[ibottom].level + 1;
				CLOSED[index] = 1;
			}
			index = bottomRules[ibottom].allSpecific.find_next(index);
		}
		//bottomRules[featureRoots_index[ibottom]].allSpecific
	}

    return bottomRules;
}


std::vector<bottomFeature> Example::initBottomFeaturesbitset(std::vector<Node *> *enrichNodes, mydataframe *enrichPropositionTable, int onlyPosFeatures)
{
    std::vector<bottomFeature> bottomRules(enrichNodes->size());
    int bitsetSize = enrichNodes->size();
    if(!onlyPosFeatures)
            bitsetSize = 2*bitsetSize;

    boost::unordered_map<Node*, int> nodePosition;
    for(int i = 0; i < bottomRules.size(); ++i)
    {
        bottomFeature tmp;        
        //Rcpp::NumericVector vctCol = (*enrichPropositionTable)[i];
        //set cover bitset
        //boost::dynamic_bitset<> coverbitset(vctCol.size(), 0);
        //for(int ibit = 0; ibit < vctCol.size(); ++ibit)
        //{
        //    if(vctCol[ibit] == 1)
        //        coverbitset[ibit] = 1;
        //}
        tmp.exampleCovered = enrichPropositionTable->data[i];//coverbitset;
        tmp.IDind = i;
        tmp.posFeature = 1;
        tmp.actScore = 0;
        tmp.pathBestScore = 0;
        tmp.noderef = (*enrichNodes)[i];
        tmp.allChilds = boost::dynamic_bitset<>(bitsetSize, 0);
        tmp.allParents = boost::dynamic_bitset<>(bitsetSize, 0);
        tmp.allSpecific = boost::dynamic_bitset<>(bitsetSize, 0);
        tmp.allGeneral = boost::dynamic_bitset<>(bitsetSize, 0);
	tmp.level = 0;

        nodePosition[tmp.noderef] = i;

        bottomRules[i] = tmp;
    }

    boost::unordered_map<Node*, int> nodePositionNEG;
    if(!onlyPosFeatures)
    {
        int origVectSize = bottomRules.size();
        bottomRules.resize(2*origVectSize);
        for(int i = 0; i < origVectSize; ++i)
        {
            bottomFeature tmp;
            //Rcpp::NumericVector vctCol = (*enrichPropositionTable)[i];
            //set cover bitset
            //boost::dynamic_bitset<> coverbitset(vctCol.size(), 0);
            //for(int ibit = 0; ibit < vctCol.size(); ++ibit)
            //{
                //negative
            //    if(vctCol[ibit] == 0)
            //        coverbitset[ibit] = 1;
            //}
            tmp.exampleCovered = ~enrichPropositionTable->data[i];//coverbitset;
            tmp.IDind = i+origVectSize;
            tmp.posFeature = 0;
            tmp.actScore = 0;
            tmp.pathBestScore = 0;
            tmp.noderef = (*enrichNodes)[i];
            tmp.allChilds = boost::dynamic_bitset<>(bitsetSize, 0);
            tmp.allParents = boost::dynamic_bitset<>(bitsetSize, 0);
            tmp.allSpecific = boost::dynamic_bitset<>(bitsetSize, 0);
            tmp.allGeneral = boost::dynamic_bitset<>(bitsetSize, 0);
	    tmp.level = 0;

            nodePositionNEG[tmp.noderef] = i+origVectSize;

            bottomRules[i+origVectSize] = tmp;
        }
    }

    //reduce unused nodes in heirarchies
    for(int iench = 0; iench < enrichNodes->size(); ++iench)
    {
        std::list<Node*> OPEN;
        OPEN.push_back((*enrichNodes)[iench]);
        Node *basedNode = OPEN.front();
        boost::unordered_map<Node*, bool> CLOSED;
        while(!OPEN.empty())
        {
            Node *actNode = OPEN.front();
            OPEN.pop_front();
            std::vector<edge> childs = actNode->relationship;
            for(int ichild = 0; ichild < childs.size(); ++ichild)
            {
                //is not a bottom and is not in CLOSED
                if(childs[ichild].end != NULL && CLOSED.find(childs[ichild].end) == CLOSED.end())
                {
                    //is the enriched node
                    if(nodePosition.find(childs[ichild].end) != nodePosition.end())
                    {
                        //relation between child and parent
                        bottomRules[nodePosition[childs[ichild].end]].allParents[nodePosition[basedNode]] = 1;
                        //relation between parent and child
                        bottomRules[nodePosition[basedNode]].allChilds[nodePosition[childs[ichild].end]] = 1;
                    }
                    else    //find another node
                    {
                        OPEN.push_back(childs[ichild].end);
                    }
                }
            }
            CLOSED[actNode] = true;
        }
    }

    if(!onlyPosFeatures)
    {
        //reduce unused nodes in heirarchies
        for(int iench = 0; iench < enrichNodes->size(); ++iench)
        {
            std::list<Node*> OPEN;
            OPEN.push_back((*enrichNodes)[iench]);
            Node *basedNode = OPEN.front();
            boost::unordered_map<Node*, bool> CLOSED;
            while(!OPEN.empty())
            {
                Node *actNode = OPEN.front();
                OPEN.pop_front();
                std::vector<edge> childs = actNode->relationship;
                for(int ichild = 0; ichild < childs.size(); ++ichild)
                {
                    //is not a bottom and is not in CLOSED
                    if(childs[ichild].end != NULL && CLOSED.find(childs[ichild].end) == CLOSED.end())
                    {
                        //is the enriched node
                        if(nodePositionNEG.find(childs[ichild].end) != nodePositionNEG.end())
                        {
                            //relation between child and parent
                            bottomRules[nodePositionNEG[childs[ichild].end]].allParents[nodePositionNEG[basedNode]] = 1;
                            //relation between parent and child
                            bottomRules[nodePositionNEG[basedNode]].allChilds[nodePositionNEG[childs[ichild].end]] = 1;
                        }
                        else    //find another node
                        {
                            OPEN.push_back(childs[ichild].end);
                        }
                    }
                }
                CLOSED[actNode] = true;
            }
        }
    }

    //complete all specific and general nodes
    for(int ibottom = 0; ibottom < bottomRules.size(); ++ibottom)
    {
        bottomRules[ibottom].allSpecific = completAllChilds(&(bottomRules[ibottom]), &bottomRules);
        bottomRules[ibottom].allGeneral = completAllParents(&(bottomRules[ibottom]), &bottomRules);
    }

//find all roots
    std::list<int> OPEN;
    for(int ibottom = 0; ibottom < bottomRules.size(); ++ibottom)
    {
		if(bottomRules[ibottom].allGeneral.none())
		{
			OPEN.push_back(ibottom);
		}
	}
	
	boost::dynamic_bitset<> CLOSED(bitsetSize, 0);
	while(!OPEN.empty())
	{
		int ibottom = OPEN.front();
		OPEN.pop_front();
		size_t index = bottomRules[ibottom].allChilds.find_first();
		
		while(index != boost::dynamic_bitset<>::npos)
		{
			if(CLOSED[index] == 0)
			{
				
				OPEN.push_back(index);
				bottomRules[index].level = bottomRules[ibottom].level + 1;
				CLOSED[index] = 1;
			}
			index = bottomRules[ibottom].allChilds.find_next(index);
		}
	}

    return bottomRules;
}

std::vector<bottomFeature> Example::initBottomFeaturesPosNeg(std::vector<Node *> *enrichNodes, std::vector<Node *> *enrichNodesNEG, Rcpp::DataFrame *enrichPropositionTable)
{
    std::vector<bottomFeature> bottomRules(enrichNodes->size() + enrichNodesNEG->size());

    boost::unordered_map<Node*, int> nodePosition;
    int bitsetSizeP = enrichNodes->size();
    int bitsetSizeN = enrichNodesNEG->size();
    for(int i = 0; i < enrichNodes->size(); ++i)
    {
        bottomFeature tmp;
        Rcpp::NumericVector vctCol = (*enrichPropositionTable)[i];
        //set cover bitset
        boost::dynamic_bitset<> coverbitset(vctCol.size(), 0);
        for(int ibit = 0; ibit < vctCol.size(); ++ibit)
        {
            if(vctCol[ibit] == 1)
                coverbitset[ibit] = 1;
        }
        tmp.exampleCovered = coverbitset;
        tmp.IDind = i;
        tmp.posFeature = 1;
        tmp.actScore = 0;
        tmp.pathBestScore = 0;
        tmp.noderef = (*enrichNodes)[i];
        tmp.allChilds = boost::dynamic_bitset<>(bottomRules.size(), 0);//bitsetSizeP, 0);
        tmp.allParents = boost::dynamic_bitset<>(bottomRules.size(), 0);//bitsetSizeP, 0);
        tmp.allSpecific = boost::dynamic_bitset<>(bottomRules.size(), 0);//bitsetSizeP, 0);
        tmp.allGeneral = boost::dynamic_bitset<>(bottomRules.size(), 0);//bitsetSizeP, 0);

        nodePosition[tmp.noderef] = i;

        bottomRules[i] = tmp;
    }

    boost::unordered_map<Node*, int> nodePositionNEG;

    int ienchNeg = 0;
    for(int i = enrichNodes->size(); i < bottomRules.size(); ++i)
    {
        //std::cout << "eichNeg: " << ienchNeg << " i : " << i << std::endl;
        bottomFeature tmp;
        Rcpp::NumericVector vctCol = (*enrichPropositionTable)[i];
        //set cover bitset
        boost::dynamic_bitset<> coverbitset(vctCol.size(), 0);
        for(int ibit = 0; ibit < vctCol.size(); ++ibit)
        {
            //negative
            if(vctCol[ibit] == 0)
                coverbitset[ibit] = 1;
        }
        tmp.exampleCovered = coverbitset;
        //std::cout << "bitset" << std::endl;
        tmp.IDind = i;
        tmp.posFeature = 0;
        tmp.actScore = 0;
        tmp.pathBestScore = 0;
        tmp.noderef = (*enrichNodesNEG)[ienchNeg];
        //std::cout << "noderef" << std::endl;
        tmp.allChilds = boost::dynamic_bitset<>(bottomRules.size(), 0);//bitsetSizeN, 0);
        tmp.allParents = boost::dynamic_bitset<>(bottomRules.size(), 0);//bitsetSizeN, 0);
        tmp.allSpecific = boost::dynamic_bitset<>(bottomRules.size(), 0);//bitsetSizeN, 0);
        tmp.allGeneral = boost::dynamic_bitset<>(bottomRules.size(), 0);//bitsetSizeN, 0);

        nodePositionNEG[tmp.noderef] = i;

        bottomRules[i] = tmp;

        ++ienchNeg;
    }

 std::cout << "sucess tmp negative nodes" << std::endl;

    //reduce unused nodes in heirarchies
    for(int iench = 0; iench < enrichNodes->size(); ++iench)
    {
        std::list<Node*> OPEN;
        OPEN.push_back((*enrichNodes)[iench]);
        Node *basedNode = OPEN.front();
        boost::unordered_map<Node*, bool> CLOSED;
        while(!OPEN.empty())
        {
            Node *actNode = OPEN.front();
            OPEN.pop_front();
            std::vector<edge> childs = actNode->relationship;
            for(int ichild = 0; ichild < childs.size(); ++ichild)
            {
                //is not a bottom and is not in CLOSED
                if(childs[ichild].end != NULL && CLOSED.find(childs[ichild].end) == CLOSED.end())
                {
                    //is the enriched node
                    if(nodePosition.find(childs[ichild].end) != nodePosition.end())
                    {
                        //relation between child and parent
                        bottomRules[nodePosition[childs[ichild].end]].allParents[nodePosition[basedNode]] = 1;
                        //relation between parent and child
                        bottomRules[nodePosition[basedNode]].allChilds[nodePosition[childs[ichild].end]] = 1;						
                    }
                    else    //find another node
                    {
                        OPEN.push_back(childs[ichild].end);
                    }
                }
            }
            CLOSED[actNode] = true;
        }
    }

std::cout << "sucess reduced positives nodes" << std::endl;

    //reduce unused nodes in heirarchies
    for(int iench = 0; iench < enrichNodesNEG->size(); ++iench)
    {
        std::list<Node*> OPEN;
        OPEN.push_back((*enrichNodesNEG)[iench]);
        Node *basedNode = OPEN.front();
        boost::unordered_map<Node*, bool> CLOSED;
        while(!OPEN.empty())
        {
            Node *actNode = OPEN.front();
            OPEN.pop_front();
            std::vector<edge> childs = actNode->relationship;
            for(int ichild = 0; ichild < childs.size(); ++ichild)
            {
                //is not a bottom and is not in CLOSED
                if(childs[ichild].end != NULL && CLOSED.find(childs[ichild].end) == CLOSED.end())
                {
                    //is the enriched node
                    if(nodePositionNEG.find(childs[ichild].end) != nodePositionNEG.end())
                    {
                        //relation between child and parent
                        bottomRules[nodePositionNEG[childs[ichild].end]].allParents[nodePositionNEG[basedNode]] = 1;
                        //relation between parent and child
                        bottomRules[nodePositionNEG[basedNode]].allChilds[nodePositionNEG[childs[ichild].end]] = 1;
                                                //std::cout << "max: " << bottomRules.size() << " base: " << nodePositionNEG[basedNode] << " max: " << bottomRules[nodePositionNEG[basedNode]].allChilds.size()  << "| child: " << nodePositionNEG[childs[ichild].end] << " max: " << bottomRules[nodePositionNEG[childs[ichild].end]].allParents.size() << std::endl;
                    }
                    else    //find another node
                    {
                        OPEN.push_back(childs[ichild].end);
                    }
                }
            }
            CLOSED[actNode] = true;
        }
    }

std::cout << "sucess reduced negative nodes" << std::endl;

    //complete all specific and general nodes
    for(int ibottom = 0; ibottom < bottomRules.size(); ++ibottom)
    {
        bottomRules[ibottom].allSpecific = completAllChilds(&(bottomRules[ibottom]), &bottomRules);
        bottomRules[ibottom].allGeneral = completAllParents(&(bottomRules[ibottom]), &bottomRules);
    }

std::cout << "sucess complete all specific and general nodes" << std::endl;

    return bottomRules;
}


boost::dynamic_bitset<> Example::completAllChilds(newComplex *feature, std::vector<bottomFeature> *bottomFeatures)
{
    if(feature->rules.size() > 0)
    {
        boost::dynamic_bitset<> childs = feature->rules[0]->allChilds;
        for(int ifeat = 1; ifeat < feature->rules.size(); ++ifeat)
        {
            //edit - not use and but OR
            childs |= feature->rules[ifeat]->allChilds;

        }
        std::list<int> OPEN;
        for(int ibit = 0; ibit < childs.size(); ++ibit)
        {
            if(childs[ibit] == 1)
                OPEN.push_back(ibit);
        }
        while(!OPEN.empty())
        {
            int nodeid = OPEN.front();
            OPEN.pop_front();
            for(int ibit = 0; ibit < (*bottomFeatures)[nodeid].allChilds.size(); ++ibit)
            {
                if((*bottomFeatures)[nodeid].allChilds[ibit] == 1)
                {
                    childs[ibit] = 1;
                    OPEN.push_back(ibit);
                }
            }
        }
        return childs;
    }
    return boost::dynamic_bitset<>(0,0);
}

boost::dynamic_bitset<> Example::completAllChilds(bottomFeature *feature, std::vector<bottomFeature> *bottomFeatures)
{

    boost::dynamic_bitset<> childs = feature->allChilds;
    std::list<int> OPEN;
    for(int ibit = 0; ibit < childs.size(); ++ibit)
    {
        if(childs[ibit] == 1)
            OPEN.push_back(ibit);
    }
    while(!OPEN.empty())
    {
        int nodeid = OPEN.front();
        OPEN.pop_front();
        for(int ibit = 0; ibit < (*bottomFeatures)[nodeid].allChilds.size(); ++ibit)
        {
            if((*bottomFeatures)[nodeid].allChilds[ibit] == 1)
            {
                childs[ibit] = 1;
                OPEN.push_back(ibit);
            }
        }
    }
    return childs;
}

boost::dynamic_bitset<> Example::completAllChildsLastRule(newComplex *feature, std::vector<bottomFeature> *bottomFeatures)
{
    if(feature->rules.size() > 0)
    {
        boost::dynamic_bitset<> childs = feature->rules.back()->allChilds;
        std::list<int> OPEN;
        for(int ibit = 1; ibit < childs.size(); ++ibit)
        {
            if(childs[ibit] == 1)
                OPEN.push_back(ibit);
        }
        while(!OPEN.empty())
        {
            int nodeid = OPEN.front();
            OPEN.pop_front();
            for(int ibit = 0; ibit < (*bottomFeatures)[nodeid].allChilds.size(); ++ibit)
            {
                if((*bottomFeatures)[nodeid].allChilds[ibit] == 1)
                {
                    childs[ibit] = 1;
                    OPEN.push_back(ibit);
                }
            }
        }
        return childs;
    }
    return boost::dynamic_bitset<>(0,0);
}

boost::dynamic_bitset<> Example::completAllParents(newComplex *feature, std::vector<bottomFeature> *bottomFeatures)
{
    if(feature->rules.size() > 0)
    {
        boost::dynamic_bitset<> parents = feature->rules[0]->allParents;
        for(int ifeat = 1; ifeat < feature->rules.size(); ++ifeat)
        {
            //edit - not use and but OR
            parents |= feature->rules[ifeat]->allParents;
        }
        std::list<int> OPEN;
        for(int ibit = 0; ibit < parents.size(); ++ibit)
        {
            if(parents[ibit] == 1)
                OPEN.push_back(ibit);
        }
        while(!OPEN.empty())
        {
            int nodeid = OPEN.front();
            OPEN.pop_front();
            for(int ibit = 0; ibit < (*bottomFeatures)[nodeid].allParents.size(); ++ibit)
            {
                if((*bottomFeatures)[nodeid].allParents[ibit] == 1)
                {
                    parents[ibit] = 1;
                    OPEN.push_back(ibit);
                }
            }
        }
        //include also generated nodes
        //for(int ifeat = 0; ifeat < feature->rules.size(); ++ifeat)
        //{
        //    parents[feature->rules[ifeat]->IDind] = 1;
        //}
        return parents;
    }
    return boost::dynamic_bitset<>(0,0);
}

boost::dynamic_bitset<> Example::completAllParents(bottomFeature *feature, std::vector<bottomFeature> *bottomFeatures)
{
    boost::dynamic_bitset<> parents = feature->allParents;
    std::list<int> OPEN;
    for(int ibit = 0; ibit < parents.size(); ++ibit)
    {
        if(parents[ibit] == 1)
            OPEN.push_back(ibit);
    }
    while(!OPEN.empty())
    {
        int nodeid = OPEN.front();
        OPEN.pop_front();
        for(int ibit = 0; ibit < (*bottomFeatures)[nodeid].allParents.size(); ++ibit)
        {
            if((*bottomFeatures)[nodeid].allParents[ibit] == 1)
            {
                parents[ibit] = 1;
                OPEN.push_back(ibit);
            }
        }
    }
    return parents;
}

boost::dynamic_bitset<> Example::getAllSpecifics(newComplex *feature)
{
    boost::dynamic_bitset<> childs = feature->rules[0]->allSpecific;
    for(int ifeat = 1; ifeat < feature->rules.size(); ++ifeat)
    {
        //edit - not use and but OR
        childs |= feature->rules[ifeat]->allSpecific;

    }
    return childs;
}

boost::dynamic_bitset<> Example::getAllGenerals(newComplex *feature)
{
    boost::dynamic_bitset<> generals = feature->rules[0]->allGeneral;
    for(int ifeat = 1; ifeat < feature->rules.size(); ++ifeat)
    {
        //edit - not use and but OR
        generals |= feature->rules[ifeat]->allGeneral;

    }
    return generals;
}

void Example::addCurrNode2CLOSED(newComplex *element2Expand, boost::dynamic_bitset<> *closed)
{
    for(int irule = 0; irule < element2Expand->rules.size(); ++irule)
    {
        (*closed)[element2Expand->rules[irule]->IDind] = 1;
    }
}

boost::dynamic_bitset<> Example::getRoots(std::vector<bottomFeature> *bottom, boost::dynamic_bitset<> *restriction)
{
    boost::dynamic_bitset<> roots(restriction->size(),0);
    for(int ibottom = 0; ibottom < bottom->size(); ++ibottom)
    {
        if((*restriction)[ibottom] == 1)
        {
            if((*bottom)[ibottom].allParents.none())
            {
                (roots)[ibottom] = 1;
            }
        }
    }
    return roots;
}

boost::dynamic_bitset<> Example::getRoots(std::vector<bottomFeature> *bottom)
{
    boost::dynamic_bitset<> roots(bottom->size(),0);
    for(int ibottom = 0; ibottom < bottom->size(); ++ibottom)
    {
        if((*bottom)[ibottom].allParents.none())
        {
            (roots)[ibottom] = 1;
        }
    }
    return roots;
}
