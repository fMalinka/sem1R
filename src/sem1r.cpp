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

    int findDescription();

    void test();

    double significantTest(int iterate, Rcpp::NumericVector tmpGene, Rcpp::NumericVector tmpSample);
    Rcpp::NumericVector getParetoScore(Rcpp::NumericVector tmpGene, Rcpp::NumericVector tmpSample);
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


    int popsize;
    int ngen;
    int kbics;
    double pcross;
    double pmut;
    double eta_c;
    double eta_m;
    double f1cutScore;
    int initPop;
    int filterTh;
    int ruleDepth;
    int onlyPosFeatures;
    double positiveTH;
    double negativeTH;
    double signTH;	//siginificance threshold
    std::string objective;
    int nrules;
    int exhaustiveTest;
    int featureSelectionMethod;
    int minLevel;

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
    std::string getPrintableRuleDetail(newComplexStat *toPrint);


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
    .method("test", &sem1R::test)
    .method("setDataset", &sem1R::setDataset)
    .method("setTestDataset", &sem1R::setTestDataset)
    .method("createROWOntology", &sem1R::createROWOntology)
    .method("createCOLOntology", &sem1R::createCOLOntology)
    .method("createTestROWOntology", &sem1R::createTestROWOntology)
    .method("createTestCOLOntology", &sem1R::createTestCOLOntology)
    .method("significantTest", &sem1R::significantTest)
    .method("getParetoScore", &sem1R::getParetoScore)
    .method("countDist", &sem1R::countDist)
    .method("combineBiclusters", &sem1R::combineBiclusters)
    .method("printColAncester", &sem1R::printColAncester)
    .method("printAllColAncester", &sem1R::printAllColAncester)
    .method("printColTopAncester", &sem1R::printColTopAncester)
    .method("getSimilaritySemanticMatrix", &sem1R::getSimilaritySemanticMatrix)
    .method("getXORscore", &sem1R::getXORscore)
    .method("getCorrectscore", &sem1R::getCorrectscore)
    .method("getSemanticRowscore", &sem1R::getSemanticRowscore)
    .method("getSemanticColscore", &sem1R::getSemanticColscore)
    .method("printParetoSettings", &sem1R::printParetoSettings)
    .method("fisher_test_cpp", &sem1R::fisher_test_cpp)
    .method("getPropositionalDataFrame", &sem1R::getPropositionalDataFrame)
    .method("getEnrichPropositionalDataFrame", &sem1R::getEnrichPropositionalDataFrame)
    .method("getEnrichPropositionalDataFrame2", &sem1R::getEnrichPropositionalDataFrame2)
    .method("getROCcordinates", &sem1R::getROCcordinates)
    .method("getCoveredVector", &sem1R::getCoveredVector)
    .method("getClassLabel", &sem1R::getClassLabel)
    .field("popsize", &sem1R::popsize)
    .field("ngen", &sem1R::ngen)
    .field("kbics", &sem1R::kbics)
    .field("pcross", &sem1R::pcross)
    .field("pmut", &sem1R::pmut)
    .field("eta_c", &sem1R::eta_c)
    .field("eta_m", &sem1R::eta_m)
    .field("f1cutScore", &sem1R::f1cutScore)
    .field("initPop", &sem1R::initPop)
    .field("filterTh", &sem1R::filterTh)
    .field("ruleDepth", &sem1R::ruleDepth)
    .field("onlyPosFeatures", &sem1R::onlyPosFeatures)
    .field("positiveTH", &sem1R::positiveTH)
    .field("negativeTH", &sem1R::negativeTH)
    .field("signTH", &sem1R::signTH)
    .field("objective", &sem1R::objective)
    .field("nrules", &sem1R::nrules)
    .field("exhaustiveTest", &sem1R::exhaustiveTest)
    .field("featureSelectionMethod", &sem1R::featureSelectionMethod)
    .field("minLevel", &sem1R::minLevel)
;
}

sem1R::sem1R()
{
    this->popsize = 200;
    this->ngen = 5000;
    this->kbics = 3;
    this->pcross = 0.8;
    this->pmut = 0.0001;
    this->eta_c = 10;
    this->eta_m = 10;
    this->f1cutScore = 0.5;
    this->initPop = 1;
    this->filterTh = 2000;
    this->ruleDepth = 2;
    this->onlyPosFeatures = 1;
    this->positiveTH = 0.7;
    this->negativeTH = 0.2;
    this->signTH = 6.635;
    this->objective = "f1";
    this->nrules = 2;
    this->exhaustiveTest = 0;
    this->featureSelectionMethod = 0;	/* 0 mincoveref, 1 hypergeometric */
    this->minLevel = -1;
}

void sem1R::printParetoSettings()
{
    Rcpp::Rcout << "[PARETO SETTINGS]" << std::endl;
    Rcpp::Rcout << "popsize: " << popsize << std::endl;
    Rcpp::Rcout << "ngen: " << ngen << std::endl;
    Rcpp::Rcout << "kbics: " << kbics << std::endl;
    Rcpp::Rcout << "pcross: " << pcross << std::endl;
    Rcpp::Rcout << "pmut: " << pmut << std::endl;
    Rcpp::Rcout << "eta_c: " << eta_c << std::endl;
    Rcpp::Rcout << "eta_m: " << eta_m << std::endl;
    Rcpp::Rcout << "f1cutScore: " << f1cutScore << std::endl;
    Rcpp::Rcout << "initPop: " << initPop << std::endl;
    Rcpp::Rcout << "filter threshold: " << filterTh << std::endl;
    Rcpp::Rcout << "rule depth: " << ruleDepth << std::endl;
    Rcpp::Rcout << "only positive features: " << onlyPosFeatures << std::endl;
    Rcpp::Rcout << "positive threshold: " << positiveTH << std::endl;
    Rcpp::Rcout << "negative threshold: " << negativeTH << std::endl;
    Rcpp::Rcout << "significance threshold: " << signTH << std::endl;
    Rcpp::Rcout << "objective function: " << objective << std::endl;
    Rcpp::Rcout << "number of rules: " << nrules << std::endl;
    Rcpp::Rcout << "exhaustive test: " << exhaustiveTest << std::endl;
    Rcpp::Rcout << "featureSelectionMethod: " << featureSelectionMethod << std::endl;
}

void sem1R::test()
{
    std::vector<Ontology *> refOntologies;

    for(boost::unordered_map<std::string, Ontology *>::iterator it = this->rowOntology.begin(); it != this->rowOntology.end(); ++it)
    {
        refOntologies.push_back(it->second);
    }
    for(boost::unordered_map<std::string, Ontology *>::iterator it = this->colOntology.begin(); it != this->colOntology.end(); ++it)
    {
        refOntologies.push_back(it->second);
    }

    //examples for the k-th bicluster
    Example bicExamples(&(this->testarmaData), &(this->rowTestOntologyDesc), &(this->colTestOntologyDesc), &(refOntologies));
    //recognize positive and negative examples
    //bicExamples.findPositiveExamplesv2(this->positiveTH);
    bicExamples.buildTestingExamples(this->positiveTH, this->negativeTH);
    //bicExamples.findNegativeExamplesv2(this->negativeTH);

//return;
    //iterate over set of rules
    boost::dynamic_bitset<> bitclass = bicExamples.getClassBitMask();
    boost::dynamic_bitset<> covered(bitclass.size(),0);
    for(int irule = 0; irule < this->ruleset.size(); ++irule)
    {

        boost::dynamic_bitset<> actCovered = bicExamples.coveredExamplesByRule(&(this->ruleset[irule]));
        covered |= actCovered;

        int tp = (covered & bitclass).count();
        int fp = (covered & ~bitclass).count();
        int tn = (~covered & ~bitclass).count();
        int fn = (~covered & bitclass).count();
        std::cout << "rule: " << irule << " tp: " << tp << " fp: " << fp << " tn: " << tn << " fn: " << fn << " acc: " << (tp+tn)/(double)(tp+tn+fp+fn) << " f1score: " << (2*tp)/(double)(2*tp +fp + fn) << std::endl << std::endl << std::endl;
    }
    for(int iii = 0; iii < covered.size(); ++iii)
    {
        std::cout << covered[iii] << ",";
    }
    std::cout << std::endl;

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


Rcpp::LogicalVector sem1R::getClassLabel()
{
    std::vector<Ontology *> refOntologies;

    for(boost::unordered_map<std::string, Ontology *>::iterator it = this->rowOntology.begin(); it != this->rowOntology.end(); ++it)
    {
        refOntologies.push_back(it->second);
    }
    for(boost::unordered_map<std::string, Ontology *>::iterator it = this->colOntology.begin(); it != this->colOntology.end(); ++it)
    {
        refOntologies.push_back(it->second);
    }

    //examples for the k-th bicluster
    Example bicExamples(&(this->testarmaData), &(this->rowTestOntologyDesc), &(this->colTestOntologyDesc), &(refOntologies));
    bicExamples.buildTestingExamples(this->positiveTH, this->negativeTH);

    //iterate over set of rules
    boost::dynamic_bitset<> bitclass = bicExamples.getClassBitMask();
    Rcpp::LogicalVector classvector(bitclass.size(), 0);

    for(int iv = 0; iv < classvector.size(); ++iv)
    {
        if(bitclass[iv] == 1)
            classvector(iv) = 1;
    }
    return classvector;
}


Rcpp::LogicalVector sem1R::getCoveredVector()
{
    std::vector<Ontology *> refOntologies;

    for(boost::unordered_map<std::string, Ontology *>::iterator it = this->rowOntology.begin(); it != this->rowOntology.end(); ++it)
    {
        refOntologies.push_back(it->second);
    }
    for(boost::unordered_map<std::string, Ontology *>::iterator it = this->colOntology.begin(); it != this->colOntology.end(); ++it)
    {
        refOntologies.push_back(it->second);
    }

    //examples for the k-th bicluster
    Example bicExamples(&(this->testarmaData), &(this->rowTestOntologyDesc), &(this->colTestOntologyDesc), &(refOntologies));
    bicExamples.buildTestingExamples(this->positiveTH, this->negativeTH);

    //iterate over set of rules
    boost::dynamic_bitset<> bitclass = bicExamples.getClassBitMask();
    boost::dynamic_bitset<> covered(bitclass.size(),0);
    for(int irule = 0; irule < this->ruleset.size(); ++irule)
    {

        boost::dynamic_bitset<> actCovered = bicExamples.coveredExamplesByRule(&(this->ruleset[irule]));
        covered |= actCovered;
    }
    Rcpp::LogicalVector rcppvector(covered.size(), 0);
    for(int ic = 0; ic < covered.size(); ++ic)
    {
        if(covered[ic] == 1)
            rcppvector(ic) = 1;
    }

    return rcppvector;
}

Rcpp::DataFrame sem1R::getROCcordinates()
{
    std::vector<Ontology *> refOntologies;

    for(boost::unordered_map<std::string, Ontology *>::iterator it = this->rowOntology.begin(); it != this->rowOntology.end(); ++it)
    {
        refOntologies.push_back(it->second);
    }
    for(boost::unordered_map<std::string, Ontology *>::iterator it = this->colOntology.begin(); it != this->colOntology.end(); ++it)
    {
        refOntologies.push_back(it->second);
    }

    //examples for the k-th bicluster
    Example bicExamples(&(this->testarmaData), &(this->rowTestOntologyDesc), &(this->colTestOntologyDesc), &(refOntologies));
    bicExamples.buildTestingExamples(this->positiveTH, this->negativeTH);

    //iterate over set of rules
    boost::dynamic_bitset<> bitclass = bicExamples.getClassBitMask();
    boost::dynamic_bitset<> covered(bitclass.size(),0);
    Rcpp::NumericVector fprate(this->ruleset.size(), 0);
    Rcpp::NumericVector tprate(this->ruleset.size(), 0);
    for(int irule = 0; irule < this->ruleset.size(); ++irule)
    {

        boost::dynamic_bitset<> actCovered = bicExamples.coveredExamplesByRule(&(this->ruleset[irule]));
        covered |= actCovered;

        int tp = (covered & bitclass).count();
        int fp = (covered & ~bitclass).count();
        int tn = (~covered & ~bitclass).count();
        int fn = (~covered & bitclass).count();

        fprate(irule) = fp/(double)(fp+tn);
        tprate(irule) = tp/(double)(tp+fn);

        //std::cout << "rule: " << irule << " tp: " << tp << " fp: " << fp << " tn: " << tn << " fn: " << fn << " acc: " << (tp+tn)/(double)(tp+tn+fp+fn) << " f1score: " << (2*tp)/(double)(2*tp +fp + fn) << std::endl << std::endl << std::endl;
    }

    return DataFrame::create(_["fpr"]= fprate, _["tpr"]= tprate);
}

/*
arma::Col<size_t> sem1R::getAssignmentsKmeans(std::vector<paretoSet> *paretoSet, int kbic)
{
    arma::mat biclusters(this->getArmaDataset()->n_rows * this->getArmaDataset()->n_cols, kbic * this->getPareto().nrow());
    arma::Col<size_t> assignments;

    //iterate over pareto set
    int ibiclusters = 0;
    for(int i = 0; i < paretoSet->size(); ++i)
    {
        //+0 - gene template
        //+1 - col template
        for(int ibic = 0; ibic < (*paretoSet)[i].chromozome.size(); ibic += 2)
        {
            //arma::mat x = (*paretoSet)[i].chromozome[ibic] * arma::trans((*paretoSet)[i].chromozome[ibic+1]);
            arma::vec concatVect = arma::vectorise((*paretoSet)[i].chromozome[ibic] * arma::trans((*paretoSet)[i].chromozome[ibic+1]));
            biclusters.col(ibiclusters) = concatVect;
            ++ibiclusters;
        }

    }
    mlpack::kmeans::KMeans<> kmeans(5000);
    kmeans.Cluster(biclusters, kbic, assignments);

    return assignments;
}
*/

Rcpp::List sem1R::fisher_test_cpp(const Rcpp::NumericMatrix& x, double conf_level = 0.95){

  // Obtain environment containing function
  Rcpp::Environment base("package:stats");

  // Make function callable from C++
  Rcpp::Function fisher_test = base["fisher.test"];

  // Call the function and receive its list output
  Rcpp::List test_out = fisher_test(Rcpp::_["x"] = x,
                                    Rcpp::_["conf.level"]  = conf_level);

  // Return test object in list structure
  return test_out;
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

     /*       //to test
            std::vector<boost::dynamic_bitset<> *> sem = newColOntology->getSemanticPatterns();
            boost::unordered_map<std::string, Node*> *ontologyMap = newColOntology->getOntologyParser()->getBottomUP();

            for(std::vector<boost::dynamic_bitset<> *>::iterator it = sem.begin(); it != sem.end(); ++it)
            {
                boost::dynamic_bitset<> *bit = (*it);
                for(int ib = 0; ib < bit->size(); ++ib)
                {
                    if((*bit)[ib] == 1)
                    {
                        //ontologyMap[]
                        for(boost::unordered_map<std::string, Node*>::iterator imap = ontologyMap->begin(); imap != ontologyMap->end(); ++imap)
                        {
                           if(imap->second != NULL && imap->second->idNum == ib)
                                Rcpp::Rcout << imap->second->id << " ";
                        }
                    }
                }
                Rcpp::Rcout << std::endl;
            }*/
    /*
            Rcpp::Rcout << "22,32:" << colOntology[name]->getSemanticDistance(22,32) << std::endl;
            Rcpp::Rcout << "42,43:" << colOntology[name]->getSemanticDistance(42,43) << std::endl;
            Rcpp::Rcout << "43,44:" << colOntology[name]->getSemanticDistance(43,44) << std::endl;
            Rcpp::Rcout << "44,45:" << colOntology[name]->getSemanticDistance(44,45) << std::endl;
            Rcpp::Rcout << "45,46:" << colOntology[name]->getSemanticDistance(45,46) << std::endl;
            Rcpp::Rcout << "46,47:" << colOntology[name]->getSemanticDistance(46,47) << std::endl;
          */  //to test end
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

/*
 * NSGAII RUN TEST
 *
 */

/*
std::vector<paretoSet> sem1R::runNSGA(int popsize, int ngen, int kbics, double pcross, double pmut, double eta_c, double eta_m, double f1cutScore, int initPop)
{
    int popsize = 200;//1200; // MUST BE DIVIDED BY 4
    double pcross_bin = 0.8;
    int bitlength = data.nrow() + data.ncol();
    double pmut_bin = 1.0/((double)popsize);//((double)bitlength);
    int ngen = 100;//13000;
    double min_binvar = 0;
    double max_binvar = bitlength;

    NSGA2 *nsga = new NSGA2(this->getArmaDataset(), popsize, ngen, kbics, pcross, pmut, eta_c, eta_m, f1cutScore, &rowOntology, &colOntology, initPop);

    std::vector<paretoSet> optimum;// = nsga->run();
    return optimum;
}
*/

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

Rcpp::NumericVector sem1R::getParetoScore(Rcpp::NumericVector tmpGene, Rcpp::NumericVector tmpSample)
{
    //fix size
    int size = Rcpp::sum(tmpGene)*Rcpp::sum(tmpSample); //bylo +

    //reference values
    int refXOR = Evaluate::evaluateXOR(&tmpGene, &tmpSample, this->getArmaDataset());
    double refSemGen = Evaluate::evaluateSemantic(&tmpGene, &(this->rowOntology));
    double refSemSam = Evaluate::evaluateSemantic(&tmpSample, &(this->colOntology));

    Rcpp::NumericVector score(4);
    score[0] = size;
    score[1] = refXOR;
    score[2] = refSemGen;
    score[3] = refSemSam;
    return score;
}

int sem1R::findDescription()
{
    //if(printAndCheckPipeline(SET_DESCRIPTION_BICS) == PIPELINE_UNACCEPTED)
    //    return PIPELINE_UNACCEPTED;

    //this->printParetoSettings();

    //assign biclusters
    //arma::Col<size_t> bicAssign = this->getAssignmentsKmeans(&(this->optimum), this->kbics);

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
    
    Rcpp::CharacterVector dataRownames;
    Rcpp::CharacterVector dataColnames;

    if(!Rf_isNull(Rcpp::rownames(this->data)) || !Rf_isNull(Rcpp::colnames(this->data)))
    {
            dataRownames = Rcpp::rownames(this->data);
            dataColnames = Rcpp::colnames(this->data);
    }
    
    for(int irule = 0; irule < krules; ++irule)
    {
		arma::uvec indicesPos = arma::find(armaDataTMP == 1);
		if(indicesPos.size() == 0)
		{
			std::cout << std::endl << "All positive example were covered!" << std::endl;
			break;
		}
        std::cout << std::endl << "INDUCING RULE #" << irule << std::endl;
        //Example bicExamples(this->kbics, &(this->armaData), &(this->optimum), &(refOntologies), &bicAssign);
        Example bicExamples(this->kbics, &(armaDataTMP), &(refOntologies));
        bicExamples.findPositiveExamples();
        bicExamples.findNegativeExamples();

        //this->propositionTable = bicExamples.getPropositionalTable(0);
        std::vector<Node *> enrichNodes;
        std::vector<Node *> enrichNodesNegatives;
        
        //preprocessing step
        //enrichment analysis
        
        
        /* OLD VERSION 
        if(this->featureSelectionMethod == 1)
        {
			this->enrichPropositionTable2 = bicExamples.getEnrichTablev2(0, &enrichNodes);
		}
		else if(this->featureSelectionMethod == 2)
		{
			this->enrichPropositionTable2 = bicExamples.getEnrichTablev2Ancestors(0, &enrichNodes, 1);
		}
		else
		{
			this->enrichPropositionTable2 = bicExamples.getEnrichTablev3TH(0, &enrichNodes, 1); //pareto bic, nodes, minth
			if(this->onlyPosFeatures == 0)
			{
				bicExamples.getEnrichTablev3THNegative(0, &enrichNodesNegatives, 1, &(this->enrichPropositionTable2));
			}
		}

        std::cout << "enrichNodes size: " << enrichNodes.size() << std::endl;
        std::cout << "enrichNodesNegative size: " << enrichNodesNegatives.size() << std::endl;
        boost::dynamic_bitset<> newclass = bicExamples.getEnrichClass(&(this->enrichPropositionTable2));
        std::vector<bottomFeature> bottomRules;
        if(this->onlyPosFeatures)
            bottomRules = bicExamples.initBottomFeatures(&enrichNodes, &(this->enrichPropositionTable2), this->onlyPosFeatures);
        else
            bottomRules = bicExamples.initBottomFeaturesPosNeg(&enrichNodes, &enrichNodesNegatives,  &(this->enrichPropositionTable2));
            * 
            */
            
            
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

        std::cout << "enrichNodes size: " << enrichNodes.size() << std::endl;
        std::cout << "enrichNodesNegative size: " << enrichNodesNegatives.size() << std::endl;
        boost::dynamic_bitset<> newclass = this->enrichPropositionTable3.predClass;//bicExamples.getEnrichClass(&(this->enrichPropositionTable2));
        std::vector<bottomFeature> bottomRules;
        //if(this->onlyPosFeatures)
            bottomRules = bicExamples.initBottomFeaturesbitset(&enrichNodes, &(this->enrichPropositionTable3), this->onlyPosFeatures);
        //else
        //    bottomRules = bicExamples.initBottomFeaturesPosNeg(&enrichNodes, &enrichNodesNegatives,  &(this->enrichPropositionTable2));


        BeamTopDown *beamTD = new BeamTopDown(&bottomRules, &newclass, this->objective);

        std::cout << std::endl << "################## TOPOLOGY VERSION ##################" << std::endl;
        clock_t begin = std::clock();
        newComplex actbest = beamTD->runTopologyVersion(this->filterTh, ruleDepth, this->signTH, this->minLevel);
		clock_t end = std::clock();  
		
		clock_t beginExh;
		clock_t endExh;
		if(this->exhaustiveTest)
		{
			BeamTopDown *beamTDExh = new BeamTopDown(&bottomRules, &newclass, this->objective);

			std::cout << std::endl << "################## EXHAUSTIVE VERSION ##################" << std::endl;
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
        std::vector<bottomFeature> tmprules;
        for(int ir = 0; ir < actbest.rules.size(); ++ir)
        {
            statActBest.rules.push_back(*(actbest.rules[ir]));
        }
        //statActBest.rules = actbest.rules; 
        statActBest.coverPos = actbest.coverPos;
        statActBest.coverNeg = actbest.coverNeg;
        this->ruleset[irule] = statActBest;

        std::string coveredDesc = beamTD->printOnlyCoveredColExamples(actbest, &(armaDataTMP), newclass, dataRownames, dataColnames);
        coveredString[irule] = coveredDesc;
					
        arma::mat new_armaData = beamTD->removeCoveredExamples(actbest, &(armaDataTMP), newclass);

        armaDataTMP = new_armaData;
        
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        std::cout << "[TOPOLOGY VERSION] elapsed overall time: " << elapsed_secs;
        if(this->exhaustiveTest)
        {
                double elapsed_secsExh = double(endExh - beginExh) / CLOCKS_PER_SEC;
                std::cout << std::endl << "[EXHAUSTIVE VERSION] elapsed overall time: " << elapsed_secsExh;
        }
    }
		
    clock_t endTOP = std::clock();
    double elapsed_secs_TOP = double(endTOP - beginTOP) / CLOCKS_PER_SEC;
    std::cout  << std::endl << "FINAL TIME TOPOLOGY:" << elapsed_secs_TOP << " k:" << krules << std::endl;
    
    std::cout << std::endl << std::endl  << "******************************************************************" << std::endl;
    std::cout << "************************** FINAL RULESET *************************" << std::endl;
    for(int irule = 0; irule < this->ruleset.size(); ++irule)
    {
		if(this->ruleset[irule].score >= 0)
		{
			std::cout << "===== RULE " << irule+1 << "=====" << std::endl;
			std::cout << " STATS: score " << this->ruleset[irule].score << " t-score: " << this->ruleset[irule].sig  << " POSITIVE: " << this->ruleset[irule].coverPos << " NEGATIVE: " << this->ruleset[irule].coverNeg << std::endl;
			std::cout << " RULE: " << getPrintableRuleID(&(this->ruleset[irule])) << std::endl;
			std::cout << " DETAILS: " << std::endl << getPrintableRuleDetail(&(this->ruleset[irule])) << std::endl;
			std::cout << "COVERED:" << std::endl;
			std::cout << coveredString[irule] << std::endl;
			std::cout << "===== =====" << std::endl << std::endl;
		}
    }




// END RULE LEARNING
//return 1;
/*
    std::cout << std::endl << std::endl << "################## EXHAUSTIVE VETSION ##################" << std::endl;
    int krules2 = 5;
    arma::mat armaDataTMP2 = this->armaData;
    clock_t beginEXH = std::clock();
    for(int irule = 0; irule < krules2; ++irule)
    {

        std::cout << std::endl << "INDUCING RULE #" << irule << std::endl;
        //Example bicExamples(this->kbics, &(this->armaData), &(this->optimum), &(refOntologies), &bicAssign);
        Example bicExamples(this->kbics, &(armaDataTMP2), &(refOntologies));
        bicExamples.findPositiveExamples();
        bicExamples.findNegativeExamples();

        //this->propositionTable = bicExamples.getPropositionalTable(0);
        std::vector<Node *> enrichNodes;
        std::vector<Node *> enrichNodesNegatives;
        //enrichment analysis
        //this->enrichPropositionTable2 = bicExamples.getEnrichTablev2(0, &enrichNodes);

        this->enrichPropositionTable2 = bicExamples.getEnrichTablev3TH(0, &enrichNodes, 1); //pareto bic, nodes, minth
        if(this->onlyPosFeatures == 0)
        {
            bicExamples.getEnrichTablev3THNegative(0, &enrichNodesNegatives, 1, &(this->enrichPropositionTable2));
        }

        std::cout << "enrichNodes size: " << enrichNodes.size() << std::endl;
        std::cout << "enrichNodesNegative size: " << enrichNodesNegatives.size() << std::endl;
        boost::dynamic_bitset<> newclass = bicExamples.getEnrichClass(&(this->enrichPropositionTable2));
        std::vector<bottomFeature> bottomRules;
        if(this->onlyPosFeatures)
            bottomRules = bicExamples.initBottomFeatures(&enrichNodes, &(this->enrichPropositionTable2), this->onlyPosFeatures);
        else
            bottomRules = bicExamples.initBottomFeaturesPosNeg(&enrichNodes, &enrichNodesNegatives,  &(this->enrichPropositionTable2));

        BeamTopDown *beamTD = new BeamTopDown(&bottomRules, &newclass);
        clock_t begin = std::clock();
        newComplex actbest = beamTD->runExhaustive(this->filterTh, ruleDepth);
        arma::mat new_armaData = beamTD->removeCoveredExamples(actbest, &(armaDataTMP2));
        armaDataTMP2 = new_armaData;

        clock_t end = std::clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        std::cout << "elapsed overall time: " << elapsed_secs;
    }

    clock_t endEXH = std::clock();
    double elapsed_secs_EXH = double(endEXH - beginEXH) / CLOCKS_PER_SEC;
    std::cout << "FINAL TIME EXHAUSTED:" << elapsed_secs_EXH << " k:" << krules2 << std::endl;






    //begin = std::clock();

    //beamTD->runExhaustive(this->filterTh, ruleDepth);//&(this->enrichPropositionTable));
    //end = std::clock();
    //elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    //std::cout << "elapsed overall time: " << elapsed_secs;
*/
    return 1;


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



void sem1R::combineBiclusters(Rcpp::NumericVector tmpGene1, Rcpp::NumericVector tmpSample1, Rcpp::NumericVector tmpGene2, Rcpp::NumericVector tmpSample2, Rcpp::NumericMatrix data)
{
 /*   Rcpp::NumericVector geneDiffer = Rcpp::abs(tmpGene1 - tmpGene2);
    Rcpp::NumericVector sampleDiffer = Rcpp::abs(tmpSample1 - tmpSample2);

    Rcpp::NumericVector geneOne(1, geneDiffer.size());
    Rcpp::NumericVector sampleOne(1, sampleDiffer.size());
    Rcpp::NumericVector geneCore = Rcpp::abs(geneDiffer - geneOne);
    Rcpp::NumericVector sampleCore = Rcpp::abs(sampleDiffer - sampleOne);

    std::vector<int> geneDifferID;
    std::vector<int> sampleDifferID;
    for(int i = 0; i < geneDiffer.size(); ++i)
    {
        if(geneDiffer(i) == 1)
            geneDifferID.push_back(i);
    }

    for(int i = 0; i < sampleDiffer.size(); ++i)
    {
        if(sampleDiffer(i) == 1)
            sampleDifferID.push_back(i);
    }

    //only genes differ
    bool sampleEqual = Rcpp::all(tmpSample1 == tmpSample2).is_true();
    bool geneEqual = Rcpp::all(tmpGene1 == tmpGene2).is_true();
    if(sampleEqual && geneEqual)
    {
        //Rcpp::Rcout << " Biclusters are equal" << std::endl;
        int size = Rcpp::sum(tmpGene1)*Rcpp::sum(tmpSample1);

        //reference values
        int refXOR = Evaluate::evaluateXOR(&tmpGene1, &tmpSample1, &data);
        double refSemGen = Evaluate::evaluateSemantic(&tmpGene1, &(this->rowOntology));
        double refSemSam = Evaluate::evaluateSemantic(&tmpSample1, &(this->colOntology));
        Rcpp::Rcout << "allequal - size: " << size << " refXOR:" << refXOR << " refSemGen:" << refSemGen << " refSemSam:" << refSemSam << std::endl;
    }
    else if(sampleEqual)
    {
        for(int igene = 0; igene < geneDifferID.size(); ++igene)
        {
            Rcpp::NumericVector tmpGene = geneCore;
            Rcpp::NumericVector tmpSample = sampleCore;
            tmpGene(geneDifferID[igene]) = 1;

            int size = Rcpp::sum(tmpGene)*Rcpp::sum(tmpSample);

            //reference values
            int refXOR = Evaluate::evaluateXOR(&tmpGene, &tmpSample, &data);
            double refSemGen = Evaluate::evaluateSemantic(&tmpGene, &(this->rowOntology));
            double refSemSam = Evaluate::evaluateSemantic(&tmpSample, &(this->colOntology));
            Rcpp::Rcout << "samplequal - size: " << size << " refXOR:" << refXOR << " refSemGen:" << refSemGen << " refSemSam:" << refSemSam << std::endl;
        }
    }
    //only samples differ
    else if(geneEqual)
    {
        for(int isample = 0; isample < sampleDifferID.size(); ++isample)
        {
            Rcpp::NumericVector tmpGene = geneCore;
            Rcpp::NumericVector tmpSample = sampleCore;
            tmpSample(sampleDifferID[isample]) = 1;

            int size = Rcpp::sum(tmpGene)*Rcpp::sum(tmpSample);

            //reference values
            int refXOR = Evaluate::evaluateXOR(&tmpGene, &tmpSample, &data);
            double refSemGen = Evaluate::evaluateSemantic(&tmpGene, &(this->rowOntology));
            double refSemSam = Evaluate::evaluateSemantic(&tmpSample, &(this->colOntology));
            Rcpp::Rcout << "genequal - size: " << size << " refXOR:" << refXOR << " refSemGen:" << refSemGen << " refSemSam:" << refSemSam << std::endl;
        }
    }
    else
    {
        for(int igene = 0; igene < geneDifferID.size(); ++igene)
        {
            for(int isample = 0; isample < sampleDifferID.size(); ++isample)
            {
                Rcpp::NumericVector tmpGene = geneCore;
                Rcpp::NumericVector tmpSample = sampleCore;
                tmpGene(geneDifferID[igene]) = 1;
                tmpSample(sampleDifferID[isample]) = 1;

                int size = Rcpp::sum(tmpGene)*Rcpp::sum(tmpSample);

                //reference values
                int refXOR = Evaluate::evaluateXOR(&tmpGene, &tmpSample, &data);
                double refSemGen = Evaluate::evaluateSemantic(&tmpGene, &(this->rowOntology));
                double refSemSam = Evaluate::evaluateSemantic(&tmpSample, &(this->colOntology));
                Rcpp::Rcout << "unequal - size: " << size << " refXOR:" << refXOR << " refSemGen:" << refSemGen << " refSemSam:" << refSemSam << std::endl;
            }
        }
    }
*/
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
