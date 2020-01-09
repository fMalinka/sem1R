//#include <Rcpp.h>

#include "RcppMLPACK.h"
#include "ontologyParser.h"
#include "ontology.h"
#include "hypothesis.h"
// #include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include "Evaluate.h"
#include <ctime>

#include "biclusterDesc/beamTopDown/conjunction_learning.h"
#include "biclusterDesc/beamTopDown/example.h"
#include "biclusterDesc/beamTopDown/ontologyrules.h"

#include <progress.hpp>
#include <progress_bar.hpp>

using namespace Rcpp;

class sem1R {
public:
    //constructors
    sem1R();

    void createROWOntology(std::string name, std::string pathToOntology, Rcpp::List descriptionTerms);
    void createCOLOntology(std::string name, std::string pathToOntology, Rcpp::List descriptionTerms);
    void createTestROWOntology(std::string name, Rcpp::List descriptionTerms);
    void createTestCOLOntology(std::string name, Rcpp::List descriptionTerms);

    void loadDatasets();

    void setDataset(Rcpp::NumericMatrix data_) {this->data = data_; armaData = as<arma::mat>(data_); pipelineCODE[SET_DATASET] = 1;}
    void setTestDataset(Rcpp::NumericMatrix data_) {this->testdata = data_; testarmaData = as<arma::mat>(data_);}
    Rcpp::NumericMatrix getDataset() { return data; }
    arma::mat* getArmaDataset() {return &(this->armaData);}

    Rcpp::List findDescription();

    double significantTest(int iterate, Rcpp::NumericVector tmpGene, Rcpp::NumericVector tmpSample);    
    Rcpp::NumericMatrix countDist(Rcpp::NumericMatrix bicMatrix);
    void combineBiclusters(Rcpp::NumericVector tmpGene1, Rcpp::NumericVector tmpSample1, Rcpp::NumericVector tmpGene2, Rcpp::NumericVector tmpSample2, Rcpp::NumericMatrix data);
    //void combineBiclusters(Rcpp::NumericVector tmpGene1, Rcpp::NumericVector tmpSample1, Rcpp::NumericVector *tmpGene2, Rcpp::NumericVector *tmpSample2, Rcpp::NumericMatrix *data);

    void printColAncester(std::string fbgnid);
    void printAllColAncester(std::string fbgnid);
    void printColTopAncester(std::string fbgnid);
    Rcpp::NumericMatrix getSimilaritySemanticMatrix(std::string ontologyName);
    int getXORscore(Rcpp::NumericVector tmpGene, Rcpp::NumericVector tmpSample);
    int getCorrectscore(Rcpp::NumericVector tmpGene, Rcpp::NumericVector tmpSample);
    double getSemanticRowscore(Rcpp::NumericVector tmp);
    double getSemanticColscore(Rcpp::NumericVector tmp);
    void printParetoSettings();

    Rcpp::List fisher_test_cpp(const Rcpp::NumericMatrix& x, double conf_level);


    Rcpp::DataFrame getPropositionalDataFrame();
    Rcpp::DataFrame getEnrichPropositionalDataFrame();
    Rcpp::DataFrame getEnrichPropositionalDataFrame2();

    Rcpp::LogicalVector getClassLabel();

    Rcpp::DataFrame getROCcordinates();
    Rcpp::LogicalVector getCoveredVector();

    int exhaustiveTest;
    int filterTh;
    int ruleDepth;
    int onlyPosFeatures;
    double signTH;	//siginificance threshold
    std::string objective;
    int nrules;    
    int featureSelectionMethod;
    int minLevel;
    std::string ruleFormat;
    int verbose;

private:
    boost::unordered_map<std::string, Ontology*> rowOntology;
    boost::unordered_map<std::string, Ontology*> colOntology;    
    boost::unordered_map<std::string, std::vector<std::vector<std::string> > > rowTestOntologyDesc;
    boost::unordered_map<std::string, std::vector<std::vector<std::string> > > colTestOntologyDesc;

    Rcpp::NumericMatrix data;
    Rcpp::NumericMatrix testdata;
    arma::mat armaData;
    arma::mat testarmaData;
    //std::vector<paretoSet> optimum;

    std::vector<newComplexStat> ruleset;

    Rcpp::DataFrame propositionTable;
    Rcpp::DataFrame enrichPropositionTable;
    Rcpp::DataFrame enrichPropositionTable2;
    mydataframe enrichPropositionTable3;

    int checkTermDescription();
    void printHypothesis(hypothesis *hyp);
    //double evaluate(arma::vec *rowBitset, arma::vec *colBitset, hypothesis *hyp, double (sem1R::*fcePtrEval)(arma::mat *, arma::mat *));
    std::vector<hypothesis> findNewLGN();
    //arma::Col<size_t> getAssignmentsKmeans(std::vector<paretoSet> *paretoSet, int kbic);
    int getPosExamples();
    int getNegExamples();
    double getHyperGeometricScore(int *successSample, int *successPop, int *failurePop, int *sampleSize);

    //evaluation functions
    double myVAR(arma::mat *templ, arma::mat *original);
    double myXOR(arma::mat *templ, arma::mat *original);
    
    std::string getPrintableRuleID(newComplexStat *toPrint);
    Rcpp::CharacterVector getRuleID(newComplexStat *toPrint);
    std::string getPrintableRuleDetail(newComplexStat *toPrint);
    Rcpp::CharacterVector getRuleDetail(newComplexStat *toPrint);


    //std::vector<paretoSet> runNSGA(int population, int generations, int kbics, double pcross, double pmut, double eta_c, double eta_m, double f1cutScore, int initPop);
    void buildTrainingExamples(std::vector< std::vector<boost::dynamic_bitset<>* > > *examples, arma::umat *indexes, std::vector<Ontology *> *refOntologies);
    void buildTrainingExamples(std::vector< std::vector<boost::dynamic_bitset<>* > > *examples, arma::umat *indexMatrix, std::vector<Ontology *> *refOntologies, std::vector<boost::unordered_map<int, bool> > *indexPosExamples);
};
RCPP_MODULE(mod_data)
{
    class_<sem1R>("sem1R")

    // expose the default constructor
    .constructor()

    .method("findDescription", &sem1R::findDescription)
    .method("setDataset", &sem1R::setDataset)
    .method("setTestDataset", &sem1R::setTestDataset)
    .method("createROWOntology", &sem1R::createROWOntology)
    .method("createCOLOntology", &sem1R::createCOLOntology)
    .method("createTestROWOntology", &sem1R::createTestROWOntology)
    .method("createTestCOLOntology", &sem1R::createTestCOLOntology)
    .method("significantTest", &sem1R::significantTest)    
    .field("filterTh", &sem1R::filterTh)
    .field("ruleDepth", &sem1R::ruleDepth)
    .field("signTH", &sem1R::signTH)
    .field("objective", &sem1R::objective)
    .field("nrules", &sem1R::nrules)
    .field("exhaustiveTest", &sem1R::exhaustiveTest)
    .field("featureSelectionMethod", &sem1R::featureSelectionMethod)
    .field("minLevel", &sem1R::minLevel)
    .field("ruleFormat", &sem1R::ruleFormat)
    .field("verbose", &sem1R::verbose)
;
}

sem1R::sem1R()
{    
    this->filterTh = 2000;
    this->ruleDepth = 2;
    this->onlyPosFeatures = 1;
    this->signTH = 6.635;
    this->objective = "f1";
    this->nrules = 2;
    this->exhaustiveTest = 0;
    this->featureSelectionMethod = 0;	/* 0 mincoveref, 1 hypergeometric */
    this->minLevel = -1;
    this->ruleFormat = "both";
    this->verbose = 0;
}

void sem1R::printParetoSettings()
{
    Rcpp::Rcout << "[sem1R SETTINGS]" << std::endl;
    Rcpp::Rcout << "filter threshold: " << filterTh << std::endl;
    Rcpp::Rcout << "rule depth: " << ruleDepth << std::endl;        
    Rcpp::Rcout << "significance threshold: " << signTH << std::endl;
    Rcpp::Rcout << "objective function: " << objective << std::endl;
    Rcpp::Rcout << "number of rules: " << nrules << std::endl;
    Rcpp::Rcout << "exhaustive test: " << exhaustiveTest << std::endl;
    Rcpp::Rcout << "featureSelectionMethod: " << featureSelectionMethod << std::endl;
    Rcpp::Rcout << "ruleFormat: " << ruleFormat << std::endl;
}


std::string sem1R::getPrintableRuleID(newComplexStat *toPrint)
{
    std::string rule = "";
    if(toPrint->rules.size() > 0)
    {
        if(toPrint->rules[0].posFeature)
            rule += toPrint->rules[0].noderef->id;
        else
            rule += "!" + toPrint->rules[0].noderef->id;
        for(int irule = 1; irule < toPrint->rules.size(); ++irule)
        {
            //negative
            if(toPrint->rules[irule].posFeature)
                rule += " AND " + toPrint->rules[irule].noderef->id;
            else
                rule += " AND !" + toPrint->rules[irule].noderef->id;
        }
    }
    return rule;
}

Rcpp::CharacterVector sem1R::getRuleID(newComplexStat *toPrint)
{
    Rcpp::CharacterVector rules(toPrint->rules.size());
    if(toPrint->rules.size() > 0)
    {
        for(int irule = 0; irule < toPrint->rules.size(); ++irule)
        {
            //negative
            if(toPrint->rules[irule].posFeature)
            {
                rules[irule] = toPrint->rules[irule].noderef->id;
            }
            else
            {
                rules[irule] = "!" + toPrint->rules[irule].noderef->id;
            }
        }
    }
    return rules;
}

std::string sem1R::getPrintableRuleDetail(newComplexStat *toPrint)
{
    std::string rule = "";
    if(toPrint->rules.size() > 0)
    {
        rule += "ID: " + toPrint->rules[0].noderef->id + "\n";
        rule += "NAME: " + toPrint->rules[0].noderef->name + "\n";
        rule += "DEF: " + toPrint->rules[0].noderef->def + "\n";
        rule += "level: " + boost::lexical_cast<string>(toPrint->rules[0].level) + "\n";
        for(int irule = 1; irule < toPrint->rules.size(); ++irule)
        {
            rule += "\n";
            rule += "ID: " + toPrint->rules[irule].noderef->id + "\n";
            rule += "NAME: " + toPrint->rules[irule].noderef->name + "\n";
            rule += "DEF: " + toPrint->rules[irule].noderef->def + "\n";
            rule += "level: " + boost::lexical_cast<string>(toPrint->rules[irule].level) + "\n";
        }
    }
    return rule;
}

Rcpp::CharacterVector sem1R::getRuleDetail(newComplexStat *toPrint)
{
    Rcpp::CharacterVector rule;
    if(toPrint->rules.size() > 0)
    {
        for(int irule = 0; irule < toPrint->rules.size(); ++irule)
        {
                        //rule += "\n";
                        rule.push_back("ID: " + toPrint->rules[irule].noderef->id);
                        rule.push_back("NAME: " + toPrint->rules[irule].noderef->name);
                        rule.push_back("DEF: " + toPrint->rules[irule].noderef->def);
        }
    }
    return rule;
}

double sem1R::getHyperGeometricScore(int *successSample, int *successPop, int *failurePop, int *sampleSize)
{
    //http://users.unimi.it/marray/2007/material/day4/Lecture7.pdf
    int log = 0;
    int lower = 0;
    return R::phyper(*successSample - 1, *successPop, *failurePop, *sampleSize, lower, log );
}

Rcpp::DataFrame sem1R::getPropositionalDataFrame()
{
    return this->propositionTable;
}

Rcpp::DataFrame sem1R::getEnrichPropositionalDataFrame()
{
    return this->enrichPropositionTable;
}

Rcpp::DataFrame sem1R::getEnrichPropositionalDataFrame2()
{
    return this->enrichPropositionTable2;
}

void sem1R::printColAncester(std::string fbgnid)
{
    boost::unordered_map<std::string, Node*> *hash = this->colOntology["FBGN"]->getOntologyParser()->getBottomUP();
    Node *nnodes = (*hash)[fbgnid];
    if(nnodes->relationship.empty())
    {
        Rcpp::Rcout << "#root#" << std::endl;
    }
    Rcpp::Rcout << "node relationship size: " << nnodes->relationship.size() << std::endl;
    for(std::vector<edge>::iterator iedge = nnodes->relationship.begin(); iedge != nnodes->relationship.end(); ++iedge)
    {
        if(iedge->end != NULL)
        {
           Rcpp::Rcout << "id:" << iedge->end->id << " id2: " << iedge->end->idNum << std::endl;
        }
        else
            Rcpp::Rcout << "konec" << std::endl;
    }
}

void sem1R::printAllColAncester(std::string fbgnid)
{
    boost::unordered_map<std::string, Node*> *hash = this->colOntology["FBGN"]->getOntologyParser()->getBottomUP();
    Node *nnodes = (*hash)[fbgnid];
    if(nnodes->relationship.empty())
    {
        ;//Rcpp::Rcout << "#root#" << std::endl;
    }
    //Rcpp::Rcout << "node relationship size: " << nnodes->relationship.size() << std::endl;
    for(std::vector<edge>::iterator iedge = nnodes->relationship.begin(); iedge != nnodes->relationship.end(); ++iedge)
    {
        if(iedge->end != NULL)
        {
           Rcpp::Rcout << fbgnid << "->" << "id:" << iedge->end->id << " id2: " << iedge->end->idNum << std::endl;
           printAllColAncester(iedge->end->id);
        }
        else
            ;//Rcpp::Rcout << "konec" << std::endl;
    }
}

void sem1R::printColTopAncester(std::string fbgnid)
{
    boost::unordered_map<std::string, Node*> *hash = this->colOntology["FBGN"]->getOntologyParser()->getTopDown();
    Node *nnodes = (*hash)[fbgnid];
    if(nnodes->relationship.empty())
    {
        Rcpp::Rcout << "#root#" << std::endl;
    }
    Rcpp::Rcout << "node relationship size: " << nnodes->relationship.size() << std::endl;
    for(std::vector<edge>::iterator iedge = nnodes->relationship.begin(); iedge != nnodes->relationship.end(); ++iedge)
    {
        if(iedge->end != NULL)
        {
           Rcpp::Rcout << "id:" << iedge->end->id << " id2: " << iedge->end->idNum << std::endl;
        }
        else
            Rcpp::Rcout << "konec" << std::endl;
    }
}

void sem1R::createCOLOntology(std::string name, std::string pathToOntology, Rcpp::List descriptionTerms)
{
    if(printAndCheckPipeline(SET_COL_ONTOLOGY) == PIPELINE_UNACCEPTED)
        return;

    Ontology *newColOntology = new Ontology(name, pathToOntology, descriptionTerms, COL_ONTOLOGY);
    if(newColOntology->isCorrect())
    {
        if(colOntology[name] == NULL)
        {
            newColOntology->precomputeSemanticPatterns(this->data.ncol());
            colOntology[name] = newColOntology;
        }
        else
            Rcpp::Rcerr << "Ontology name is not unique!" << std::endl;
    }
    else
         Rcpp::Rcerr << "Ontology cannot be loaded!" << std::endl;
}

void sem1R::createROWOntology(std::string name, std::string pathToOntology, Rcpp::List descriptionTerms)
{
    if(printAndCheckPipeline(SET_ROW_ONTOLOGY) == PIPELINE_UNACCEPTED)
        return;

    Ontology *newRowOntology = new Ontology(name, pathToOntology, descriptionTerms, ROW_ONTOLOGY);
    if(newRowOntology->isCorrect())
    {
        if(rowOntology[name] == NULL)
        {
            newRowOntology->precomputeSemanticPatterns(this->data.nrow());
            rowOntology[name] = newRowOntology;
        }
        else
            Rcpp::Rcerr << "Ontology name is not unique!" << std::endl;
    }
    else
        Rcpp::Rcerr << "Ontology cannot be loaded!" << std::endl;
}

void sem1R::createTestROWOntology(std::string name, Rcpp::List descriptionTerms)
{
    //if(printAndCheckPipeline(SET_ROW_ONTOLOGY) == PIPELINE_UNACCEPTED)
    //    return;

    //Ontology *newRowOntology = new Ontology(name, pathToOntology, descriptionTerms, ROW_ONTOLOGY);

    if(rowOntology[name] != NULL)
    {
        //rowTestOntologyDesc[name] = rowOntology[name]->convertRList2Vector(descriptionTerms);
        rowOntology[name]->addDescriptionTermsTest(rowOntology[name]->convertRList2Vector(descriptionTerms));
        rowOntology[name]->precomputeSemanticPatternsTest(this->testdata.nrow());
    }
    else
        Rcpp::Rcerr << "Ontology name is not found!" << std::endl;
}

void sem1R::createTestCOLOntology(std::string name, Rcpp::List descriptionTerms)
{
    //if(printAndCheckPipeline(SET_ROW_ONTOLOGY) == PIPELINE_UNACCEPTED)
    //    return;

    //Ontology *newRowOntology = new Ontology(name, pathToOntology, descriptionTerms, ROW_ONTOLOGY);

    if(colOntology[name] != NULL)
    {
        //colTestOntologyDesc[name] = colOntology[name]->convertRList2Vector(descriptionTerms);
        colOntology[name]->addDescriptionTermsTest(colOntology[name]->convertRList2Vector(descriptionTerms));
        colOntology[name]->precomputeSemanticPatternsTest(this->testdata.ncol());
    }
    else
        Rcpp::Rcerr << "Ontology name is not found!" << std::endl;
}

Rcpp::NumericMatrix sem1R::getSimilaritySemanticMatrix(std::string ontologyName)
{
    Ontology* onto;
    if(this->rowOntology.count(ontologyName))
    {
        onto = this->rowOntology[ontologyName];
    }
    else if(this->colOntology.count(ontologyName))
    {
        onto = this->colOntology[ontologyName];
    }
    else
    {
        Rcpp::Rcerr << "Ontology name is not found!" << std::endl;
        Rcpp::NumericMatrix simMatrix;
        return simMatrix;
    }

    std::vector<std::vector<double> > matrix = onto->getSemanticSimilarityMatrix();
    Rcpp::NumericMatrix simMatrix(matrix.size(), matrix.size());
    for(int i = 0; i < matrix.size(); ++i)
    {
        for(int j = 0; j < matrix.size(); ++j)
        {
            simMatrix(i,j) = matrix[i][j];
            //Rcpp::Rcout << matrix[i][j]<< std::endl;
        }
    }
    return simMatrix;
}

int sem1R::getXORscore(Rcpp::NumericVector tmpGene, Rcpp::NumericVector tmpSample)
{
     return Evaluate::evaluateXOR(&tmpGene, &tmpSample, this->getArmaDataset());
}

int sem1R::getCorrectscore(Rcpp::NumericVector tmpGene, Rcpp::NumericVector tmpSample)
{
     return Evaluate::evaluateCorrectPrediction(&tmpGene, &tmpSample, this->getArmaDataset());
}

double sem1R::getSemanticRowscore(Rcpp::NumericVector tmp)
{
     return Evaluate::evaluateSemantic(&tmp, &(this->rowOntology));
}

double sem1R::getSemanticColscore(Rcpp::NumericVector tmp)
{
     return Evaluate::evaluateSemantic(&tmp, &(this->colOntology));
}

double sem1R::significantTest(int iterate, Rcpp::NumericVector tmpGene, Rcpp::NumericVector tmpSample)
{
    //fix size
    int size = Rcpp::sum(tmpGene)+Rcpp::sum(tmpSample);
    //set seed
    Rcpp::Environment base_env("package:base");
    Rcpp::Function set_seed_r = base_env["set.seed"];
    set_seed_r(123);

    //reference values
    int refXOR = Evaluate::evaluateXOR(&tmpGene, &tmpSample, this->getArmaDataset());
    int refCorrectPred = Evaluate::evaluateCorrectPrediction(&tmpGene, &tmpSample, this->getArmaDataset());
    double refSemGen = Evaluate::evaluateSemantic(&tmpGene, &(this->rowOntology));
    double refSemSam = Evaluate::evaluateSemantic(&tmpSample, &(this->colOntology));

    Rcpp::Rcout << "ref XOR " << refXOR << std::endl;
    Rcpp::Rcout << "ref CorrectPred " << refCorrectPred << std::endl;
    Rcpp::Rcout << "ref Sem Gen " << refSemGen << std::endl;
    Rcpp::Rcout << "ref Sem Sam " << refSemSam << std::endl;



    //Rcpp::Rcout << tmpGene.length() << std::endl;
    //Rcpp::Rcout << tmpSample.length() << std::endl;

    int betterRandom = 0;
    Rcpp::IntegerVector vectorLimit = Rcpp::seq_len(tmpGene.length()+tmpSample.length());
    for(int i = 0; i < iterate; ++i)
    {
        Rcpp::NumericVector randomBic = RcppArmadillo::sample((Rcpp::NumericVector)vectorLimit, size, FALSE); //interval, how many, replacement
        //Rcpp::Rcout << randomBic << std::endl;

        Rcpp::NumericVector rndGene(tmpGene.length(),0);
        Rcpp::NumericVector rndSample(tmpSample.length(),0);
        int nrow = rndGene.length();
        int ncol = rndSample.length();

        //make a subsets
        //iterate over rows/genes
        for(int bitrow = 0; bitrow < nrow; ++bitrow)
        {
            //for more speed use [] instead of () -> bound check will be disable
            if(randomBic[bitrow] == 1)
                rndGene(bitrow) = 1;
        }
        //iterate over columns/samples
        for(int bitcol = nrow; bitcol < nrow + ncol; ++bitcol)
        {
            //for more speed use [] instead of () -> bound check will be disable
            if(randomBic[bitcol] == 1)
                rndSample(bitcol-nrow) = 1;
        }

        if(Evaluate::evaluateCorrectPrediction(&rndGene, &rndSample, this->getArmaDataset()) <=  refCorrectPred &&
                Evaluate::evaluateSemantic(&rndGene, &(this->rowOntology)) <= refSemGen &&
                Evaluate::evaluateSemantic(&rndSample, &(this->colOntology)) <= refSemSam)
            betterRandom++;

    }
    return (double)betterRandom/iterate;
}

Rcpp::List sem1R::findDescription()
{
    //if(printAndCheckPipeline(SET_DESCRIPTION_BICS) == PIPELINE_UNACCEPTED)
    //    return PIPELINE_UNACCEPTED;

    //this->printParetoSettings();
    //prepare ontologies
    std::vector<Ontology *> refOntologies;

    for(boost::unordered_map<std::string, Ontology *>::iterator it = this->rowOntology.begin(); it != this->rowOntology.end(); ++it)
    {
        refOntologies.push_back(it->second);
    }
    for(boost::unordered_map<std::string, Ontology *>::iterator it = this->colOntology.begin(); it != this->colOntology.end(); ++it)
    {
        refOntologies.push_back(it->second);
    }

    //separate biclusters to positive examples

//RULE LEARNING

    int krules = this->nrules;
    newComplexStat initStats;
    initStats.score = -1;
    this->ruleset.resize(krules, initStats);
    arma::mat armaDataTMP = this->armaData;
    clock_t beginTOP = std::clock();
    std::vector<std::string> coveredString(this->nrules);
    std::vector<Rcpp::CharacterVector> coveredRules(this->nrules);
    
    Rcpp::CharacterVector dataRownames;
    Rcpp::CharacterVector dataColnames;

    if(!Rf_isNull(Rcpp::rownames(this->data)) || !Rf_isNull(Rcpp::colnames(this->data)))
    {
            dataRownames = Rcpp::rownames(this->data);
            dataColnames = Rcpp::colnames(this->data);
    }

    printParetoSettings();

    int trueHypothesisiSize = 0;
    //the main loop
    bool display_progress = (this->verbose > 0) ? false : true;
    Progress progressBar(krules*8, display_progress);
    int progressInc = 0;
    for(int irule = 0; irule < krules; ++irule)
    {
        arma::uvec indicesPos = arma::find(armaDataTMP == 1);
        if(indicesPos.size() == 0)
        {
            if(this->verbose)
                Rcpp::Rcout << std::endl << "All positive example were covered!" << std::endl;
            break;
        }
        if(this->verbose)
            Rcpp::Rcout << std::endl << "INDUCING RULE #" << irule << std::endl;

        //step 0

        //Example bicExamples(this->kbics, &(this->armaData), &(this->optimum), &(refOntologies), &bicAssign);
        Example bicExamples(0, &(armaDataTMP), &(refOntologies), this->verbose);
        bicExamples.findPositiveExamples();

        //step 1
        if (Progress::check_abort())
            break;
        progressBar.increment();
        ++progressInc;


        bicExamples.findNegativeExamples();

        //step 2
        if (Progress::check_abort())
            break;
        progressBar.increment();
        ++progressInc;

        std::vector<Node *> enrichNodes;        
        //preprocessing step
        //enrichment analysis            
        if(this->featureSelectionMethod == 1)
        {
            this->enrichPropositionTable3 = bicExamples.getEnrichTablev2Sigbitset(0, &enrichNodes);
        }
        else if(this->featureSelectionMethod == 2)
        {
            this->enrichPropositionTable3 = bicExamples.getEnrichTablev2AncestorsSigbitset(0, &enrichNodes, 1);
        }
        else
        {
            this->enrichPropositionTable3 = bicExamples.getEnrichTablev3THbitset(0, &enrichNodes, 1); //pareto bic, nodes, minth
        }

        //step 3
        if (Progress::check_abort())
            break;
        progressBar.increment();
        ++progressInc;

        if(this->verbose)
            Rcpp::Rcout << "enrichNodes size: " << enrichNodes.size() << std::endl;
        boost::dynamic_bitset<> newclass = this->enrichPropositionTable3.predClass;//bicExamples.getEnrichClass(&(this->enrichPropositionTable2));

        //step 4
        if (Progress::check_abort())
            break;
        progressBar.increment();
        ++progressInc;


        std::vector<bottomFeature> bottomRules;
        //if(this->onlyPosFeatures)
        bottomRules = bicExamples.initBottomFeaturesbitset(&enrichNodes, &(this->enrichPropositionTable3), this->onlyPosFeatures);
        //else
        //    bottomRules = bicExamples.initBottomFeaturesPosNeg(&enrichNodes, &enrichNodesNegatives,  &(this->enrichPropositionTable2));

        //step 5
        if (Progress::check_abort())
            break;
        progressBar.increment();
        ++progressInc;

        BeamTopDown *beamTD = new BeamTopDown(&bottomRules, &newclass, this->objective, this->verbose);

        if(this->verbose)
            Rcpp::Rcout << std::endl << "################## TOPOLOGY VERSION ##################" << std::endl;
        clock_t begin = std::clock();
        newComplex actbest = beamTD->runTopologyVersion(this->filterTh, ruleDepth, this->signTH, this->minLevel);
        clock_t end = std::clock();

        //step 6
        if (Progress::check_abort())
            break;
        progressBar.increment();
        ++progressInc;

        clock_t beginExh;
        clock_t endExh;
        if(this->exhaustiveTest)
        {
                BeamTopDown *beamTDExh = new BeamTopDown(&bottomRules, &newclass, this->objective, this->verbose);
                if(this->verbose)
                    Rcpp::Rcout << std::endl << "################## EXHAUSTIVE VERSION ##################" << std::endl;
                beginExh = std::clock();
                newComplex actbestExh = beamTDExh->runExhaustive(this->filterTh, ruleDepth, this->signTH);
                endExh = std::clock();
        }
		      
        //nenasel jsem zadne pravidlo splnujici pozadavky
        if(actbest.score < 0)
            break;
        
        newComplexStat statActBest;
        statActBest.pathScore = actbest.pathScore;
        statActBest.score = actbest.score;
        statActBest.sig = actbest.sig;

        for(int ir = 0; ir < actbest.rules.size(); ++ir)
        {
            statActBest.rules.push_back(*(actbest.rules[ir]));
        }
        //statActBest.rules = actbest.rules; 
        statActBest.coverPos = actbest.coverPos;
        statActBest.coverNeg = actbest.coverNeg;
        this->ruleset[irule] = statActBest;

        if(this->ruleset[irule].score >= 0)
        {
            ++trueHypothesisiSize;
        }

        std::string coveredDesc = beamTD->printOnlyCoveredColExamples(actbest, &(armaDataTMP), newclass, dataRownames, dataColnames);
        Rcpp::CharacterVector coveredList;
        if(this->ruleFormat == "col")
            coveredList = beamTD->getOnlyCoveredColExamples(actbest, &(armaDataTMP), dataColnames);
        else if(this->ruleFormat == "row")
            coveredList = beamTD->getOnlyCoveredRowExamples(actbest, &(armaDataTMP), dataRownames);
        else if(this->ruleFormat == "both")
            coveredList = beamTD->getOnlyCoveredRowColExamples(actbest, &(armaDataTMP), dataRownames, dataColnames);
        else
            ;

        coveredRules[irule] = coveredList;
        coveredString[irule] = coveredDesc;
					
        arma::mat new_armaData = beamTD->removeCoveredExamples(actbest, &(armaDataTMP), newclass);

        armaDataTMP = new_armaData;
        
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        //Rcpp::Rcout << "[TOPOLOGY VERSION] elapsed overall time: " << elapsed_secs;
        if(this->exhaustiveTest)
        {
                double elapsed_secsExh = double(endExh - beginExh) / CLOCKS_PER_SEC;
                //Rcpp::Rcout << std::endl << "[EXHAUSTIVE VERSION] elapsed overall time: " << elapsed_secsExh;
        }
        if(statActBest.coverPos < 1)
            break;

        //step 7

        if (Progress::check_abort())
            break;
        progressBar.increment();
        ++progressInc;
    }

    //end the progressbar
    for(int iprogress = krules*8; iprogress >= progressInc; --iprogress)
    {
        progressBar.increment();
    }


    clock_t endTOP = std::clock();
    double elapsed_secs_TOP = double(endTOP - beginTOP) / CLOCKS_PER_SEC;

    if(this->verbose)
        Rcpp::Rcout  << std::endl << "FINAL TIME TOPOLOGY:" << elapsed_secs_TOP << " k:" << krules << std::endl;
    
    Rcpp::Rcout << std::endl << std::endl  << "******************************************************************" << std::endl;
    Rcpp::Rcout << "************************** FINAL RULESET *************************" << std::endl;
    Rcpp::List hypothesis(trueHypothesisiSize);
    int ihypo = 0;
    for(int irule = 0; irule < this->ruleset.size(); ++irule)
    {
		if(this->ruleset[irule].score >= 0)
		{
                        Rcpp::Rcout << "===== RULE " << irule+1 << "=====" << std::endl;
                        Rcpp::Rcout << " STATS: score " << this->ruleset[irule].score << " t-score: " << this->ruleset[irule].sig  << " POSITIVE: " << this->ruleset[irule].coverPos << " NEGATIVE: " << this->ruleset[irule].coverNeg << std::endl;
                        Rcpp::Rcout << " RULE: " << getPrintableRuleID(&(this->ruleset[irule])) << std::endl;
                        Rcpp::Rcout << " DETAILS: " << std::endl << getPrintableRuleDetail(&(this->ruleset[irule])) << std::endl;
                        Rcpp::Rcout << "COVERED:" << std::endl;
                        Rcpp::Rcout << coveredString[irule] << std::endl;
                        Rcpp::Rcout << "===== =====" << std::endl << std::endl;
		}
                //Rcpp List
                Rcpp::List listRule = Rcpp::List::create(Rcpp::_["ruleID"] = irule+1, Rcpp::_["score"] = this->ruleset[irule].score, Rcpp::_["tscore"] = this->ruleset[irule].sig,
                Rcpp::_["positiveCovered"] = this->ruleset[irule].coverPos, Rcpp::_["negativeCovered"] = this->ruleset[irule].coverNeg,
                Rcpp::_["rules"] = getRuleID(&(this->ruleset[irule])), Rcpp::_["details"] = getRuleDetail(&(this->ruleset[irule])), Rcpp::_["covered"] = coveredRules[irule]);
                hypothesis[ihypo] = listRule;
                ++ihypo;
    }

return hypothesis;

}

void sem1R::buildTrainingExamples(std::vector< std::vector<boost::dynamic_bitset<>* > > *examples, arma::umat *indexMatrix, std::vector<Ontology *> *refOntologies)
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

void sem1R::buildTrainingExamples(std::vector< std::vector<boost::dynamic_bitset<>* > > *examples, arma::umat *indexMatrix, std::vector<Ontology *> *refOntologies, std::vector<boost::unordered_map<int, bool> > *indexPosExamples)
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


Rcpp::NumericMatrix sem1R::countDist(Rcpp::NumericMatrix bicMatrix)
{
    int nbics = bicMatrix.nrow();
    int sizeBic = bicMatrix.ncol();
    Rcpp::NumericMatrix dist(nbics, nbics);
    Rcpp::Rcout << "nbics: " << nbics << std::endl;
    Rcpp::Rcout << "sizeBic: " << sizeBic << std::endl;

    //convert vector to binary set
    std::vector< boost::dynamic_bitset<> > bitsets;
    for(int i = 0; i < nbics; ++i)
    {
        boost::dynamic_bitset<> tmpbit(sizeBic, 0);
        for(int ibit = 0; ibit < sizeBic; ++ibit)
        {
            if(bicMatrix(i,ibit) == 1)
                tmpbit[ibit] = 1;
        }
        bitsets.push_back(tmpbit);
    }

    for(int i = 0; i < nbics; ++i)
    {
        for(int j = i+1; j < nbics; ++j)
        {
            boost::dynamic_bitset<> xorset;
            xorset = bitsets[i] ^ bitsets[j];
            double sim = (double) (sizeBic - xorset.count()) / sizeBic;
            dist(i,j) = sim;
            dist(j,i) = sim;
        }
    }

    return dist;
}

double sem1R::myXOR(arma::mat *templ, arma::mat *original)
{
    arma::uvec indeces = find(*templ); //find non zero --- indexes
    arma::vec elementsTempl = (*templ)(indeces);

    arma::vec origVector = (arma::vectorise(*original));
    arma::vec elementsOrig = origVector(indeces);  //vector with selected elements

    return arma::accu(arma::abs(elementsOrig-elementsTempl));
}

double sem1R::myVAR(arma::mat *templ, arma::mat *original)
{
    arma::uvec indeces = find(*templ); //find non zero --- indexes
    arma::vec origVector = (arma::vectorise(*original));
    arma::vec elements = origVector(indeces);  //vector with selected elements
    return arma::var(elements);
}


std::vector<hypothesis> sem1R::findNewLGN()
{
    std::vector<hypothesis> newRoots;

    //find new roots for row ontologies
    for(boost::unordered_map<std::string, Ontology*>::iterator irOnto = this->rowOntology.begin(); irOnto != this->rowOntology.end(); ++irOnto)
    {
        //Rcpp::Rcout << "Starting find new roots in: " << irOnto->first << std::endl;
        std::vector<Node*> tmpRoot = irOnto->second->findLGNroots();
        for(std::vector<Node*>::iterator it = tmpRoot.begin(); it != tmpRoot.end(); ++it)
        {
            hypothesis hypothes;
            hypothes.score = 0;
            hypothes.rowOntology.push_back(irOnto->second);
            boost::unordered_map<std::string, bool> id;
            id[(*it)->id] = true;
            hypothes.ids[irOnto->second] = id;
            newRoots.push_back(hypothes);
        }
    }
    //find new roots for col ontologies
    for(boost::unordered_map<std::string, Ontology*>::iterator icOnto = this->colOntology.begin(); icOnto != this->colOntology.end(); ++icOnto)
    {
        //Rcpp::Rcout << "Starting find new roots in: " << icOnto->first << std::endl;
        std::vector<Node*> tmpRoot = icOnto->second->findLGNroots();
        for(std::vector<Node*>::iterator it = tmpRoot.begin(); it != tmpRoot.end(); ++it)
        {
            hypothesis hypothes;
            hypothes.score = 0;
            hypothes.colOntology.push_back(icOnto->second);
            boost::unordered_map<std::string, bool> id;
            id[(*it)->id] = true;
            hypothes.ids[icOnto->second] = id;
            newRoots.push_back(hypothes);
        }
    }
    return newRoots;
}

void sem1R::printHypothesis(hypothesis *hyp)
{
    for(boost::unordered_map<Ontology*,boost::unordered_map<std::string, bool> >::iterator iOnto = hyp->ids.begin(); iOnto != hyp->ids.end(); ++iOnto)
    {
        for(boost::unordered_map<std::string, bool>::iterator istring = iOnto->second.begin(); istring != iOnto->second.end(); ++istring)
        {
            Rcpp:Rcout << istring->first << " ";
        }
    }
}



int sem1R::checkTermDescription()
{
    //col ontology
    bool error = false;
    for(boost::unordered_map<std::string, Ontology*>::iterator icol = colOntology.begin(); icol != colOntology.end(); ++icol)
    {
        if(icol->second->getDescriptionTerms().size() != data.ncol())
        {
            Rcpp::Rcerr << icol->first << " ontology and term description have unequal size!" << std::endl;
            error = true;
        }
    }
    //row ontology
    for(boost::unordered_map<std::string, Ontology*>::iterator irow = rowOntology.begin(); irow != rowOntology.end(); ++irow)
    {
        if(irow->second->getDescriptionTerms().size() != data.nrow())
        {
            Rcpp::Rcerr << irow->first << " ontology and term description have unequal size!" << std::endl;
            error = true;
        }
    }
    return error;
}
