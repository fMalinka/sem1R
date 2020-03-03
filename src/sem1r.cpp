//#include <Rcpp.h>

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

    void setDataset(Rcpp::NumericMatrix data_) {this->data = data_; armaData = as<arma::mat>(data_); this->rowOntology.clear(); this->colOntology.clear(); printAndCheckPipeline(SET_DATASET);}
    Rcpp::List findDescription();
    Rcpp::List computeTermsEnrichment();
    Rcpp::List checkRowDescription();
    Rcpp::List checkColDescription();
    void printParetoSettings();

    Rcpp::List getIgraphNodes() {return this->igraphNode;}
    Rcpp::List getIgraphEdges() {return this->igraphEdge;}

    int filterTh;
    int ruleDepth;
    double signTH;	//siginificance threshold
    std::string objective;
    int nrules;    
    int featureSelectionMethod;
    int minLevel;
    std::string ruleFormat;
    int verbose;

    Rcpp::List igraphNode;
    Rcpp::List igraphEdge;

private:
    boost::unordered_map<std::string, Ontology*> rowOntology;
    boost::unordered_map<std::string, Ontology*> colOntology;        

    Rcpp::NumericMatrix data;
    arma::mat armaData;
    arma::mat testarmaData;
    std::vector<newComplexStat> ruleset;

    Rcpp::DataFrame propositionTable;
    Rcpp::DataFrame enrichPropositionTable;
    Rcpp::DataFrame enrichPropositionTable2;
    mydataframe enrichPropositionTable3;
    
    std::string getPrintableRuleID(newComplexStat *toPrint);
    Rcpp::CharacterVector getRuleID(newComplexStat *toPrint);
    std::string getPrintableRuleDetail(newComplexStat *toPrint);
    Rcpp::CharacterVector getRuleDetail(newComplexStat *toPrint);   
};

RCPP_MODULE(mod_data)
{
    class_<sem1R>("sem1R")

    // expose the default constructor
    .constructor()

    .method("findDescription", &sem1R::findDescription, "Run the algorithm, find a hypothesis with given parameters.")
    .method("setDataset", &sem1R::setDataset,  "Set a two-dimensional binary matrix.")
    .method("createROWOntology", &sem1R::createROWOntology, "Append a row ontology.")
    .method("createCOLOntology", &sem1R::createCOLOntology, "Append a column ontology.")
    .method("getIgraphNodes", &sem1R::getIgraphNodes, "Returns igraph data.frame nodes.")
    .method("getIgraphEdges", &sem1R::getIgraphEdges, "Returns igraph data.frame edges.")
    .method("computeTermsEnrichment", &sem1R::computeTermsEnrichment, "Compute enrichment score for each ontological term.")
    .method("checkRowDescription", &sem1R::checkRowDescription, "Check row ontology description terms.")
    .method("checkColDescription", &sem1R::checkColDescription, "Check col ontology description terms.")
    .field("filterTh", &sem1R::filterTh, "")
    .field("ruleDepth", &sem1R::ruleDepth, "Maximal length of induced rule.")
    .field("signTH", &sem1R::signTH, "Minimal value of Chi Square.")
    .field("objective", &sem1R::objective, "")
    .field("nrules", &sem1R::nrules, "Maximum number of induced rules.")
    .field("featureSelectionMethod", &sem1R::featureSelectionMethod, "Type of Feature Selection method.")
    .field("minLevel", &sem1R::minLevel, "Minimum level of generalization that each term in a rule has to sattisfy.")
    .field("ruleFormat", &sem1R::ruleFormat, "Format of covered examples.")
    .field("verbose", &sem1R::verbose, "Print some debug information.")
;
}

sem1R::sem1R()
{    
    this->filterTh = 2000;
    this->ruleDepth = 2;
    this->signTH = 6.635;
    this->objective = "f1";
    this->nrules = 2;
    this->featureSelectionMethod = 0;	/* 0 mincoveref, 1 hypergeometric */
    this->minLevel = -1;
    this->ruleFormat = "both";
    this->verbose = 0;

    omp_set_num_threads(1);
}

void sem1R::printParetoSettings()
{
    Rcpp::Rcout << "[sem1R SETTINGS]" << std::endl;
    Rcpp::Rcout << "filter threshold: " << filterTh << std::endl;
    Rcpp::Rcout << "rule depth: " << ruleDepth << std::endl;        
    Rcpp::Rcout << "significance threshold: " << signTH << std::endl;
    Rcpp::Rcout << "objective function: " << objective << std::endl;
    Rcpp::Rcout << "number of rules: " << nrules << std::endl;
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
            if(newColOntology->isCorrect())
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
            if(newRowOntology->isCorrect())
                rowOntology[name] = newRowOntology;
        }
        else
            Rcpp::Rcerr << "Ontology name is not unique!" << std::endl;
    }
    else
        Rcpp::Rcerr << "Ontology cannot be loaded!" << std::endl;
}

Rcpp::List sem1R::computeTermsEnrichment()
{
    if(printAndCheckPipeline(SET_DESCRIPTION_BICS) == PIPELINE_UNACCEPTED)
        return Rcpp::List();

    std::vector<Ontology *> refOntologies;

    for(boost::unordered_map<std::string, Ontology *>::iterator it = this->rowOntology.begin(); it != this->rowOntology.end(); ++it)
    {
        refOntologies.push_back(it->second);
    }
    for(boost::unordered_map<std::string, Ontology *>::iterator it = this->colOntology.begin(); it != this->colOntology.end(); ++it)
    {
        refOntologies.push_back(it->second);
    }

    //RULE LEARNING
    arma::mat armaDataTMP = this->armaData;

    Rcpp::CharacterVector dataRownames;
    Rcpp::CharacterVector dataColnames;

    if(!Rf_isNull(Rcpp::rownames(this->data)) || !Rf_isNull(Rcpp::colnames(this->data)))
    {
            dataRownames = Rcpp::rownames(this->data);
            dataColnames = Rcpp::colnames(this->data);
    }
    else
    {
        if(Rf_isNull(Rcpp::rownames(this->data)))
        {
            Rcpp::Rcerr << "Rownames of the input dataset are not set!" << std::endl;
        }
        if(Rf_isNull(Rcpp::colnames(this->data)))
        {
            Rcpp::Rcerr << "Rownames of the input dataset are not set!" << std::endl;
        }
        return Rcpp::List();
    }

    printParetoSettings();

    Example bicExamples(0, &(armaDataTMP), &(refOntologies), this->verbose);
    bicExamples.findPositiveExamples();
    bicExamples.findNegativeExamples();
    return bicExamples.getEnrichTerms();
}


Rcpp::List sem1R::checkRowDescription()
{
    if(printAndCheckPipeline(SET_DESCRIPTION_BICS) == PIPELINE_UNACCEPTED)
        return Rcpp::List();

    std::vector<Ontology *> refOntologies;

    Rcpp::List ontolist(this->rowOntology.size());
    int ionto = 0;
    for(boost::unordered_map<std::string, Ontology *>::iterator it = this->rowOntology.begin(); it != this->rowOntology.end(); ++it)
    {
        //refOntologies.push_back(it->second);

        Ontology *actOnto = it->second;        
        std::vector<std::vector<std::string> > ontodescterms = actOnto->getDescriptionTerms();
        Rcpp::List listdesc(ontodescterms.size());
        for(int idesc = 0; idesc < ontodescterms.size(); ++idesc)
        {
            Rcpp::CharacterVector rowdesc;
            for(int isubdesc = 0; isubdesc < ontodescterms[idesc].size(); ++isubdesc)
            {
                //if(actOnto->getOntologyParser()->getNodesBottomUP(&(ontodescterms[idesc][isubdesc])) != NULL)
                //Rcpp::Rcout << "#"  << ontodescterms[idesc][isubdesc] << "#- ";
                if(actOnto->getOntologyParser()->getNodesBottomUP(&(ontodescterms[idesc][isubdesc])) != NULL)
                {                    
                    //Rcpp::Rcout << "(" << ontodescterms[idesc][isubdesc] << " level: " << actOnto->getLevelOfSpecialization(actOnto->getOntologyParser()->getNodesBottomUP(&(ontodescterms[idesc][isubdesc]))) << ") ";
                    rowdesc.push_back(ontodescterms[idesc][isubdesc] + std::string(", specificity level: ") + boost::lexical_cast<string>(actOnto->getLevelOfSpecialization(actOnto->getOntologyParser()->getNodesBottomUP(&(ontodescterms[idesc][isubdesc])))));
                }
                else
                {
                    //Rcpp::Rcout << "ontodesc: " << ontodescterms[idesc][isubdesc] << "not found!";
                    rowdesc.push_back(ontodescterms[idesc][isubdesc] + std::string(" not found!"));
                }
            }
            listdesc[idesc] = rowdesc;
            Rcpp::Rcout << std::endl;
        }
        ontolist[ionto] = listdesc;
        ++ionto;
    }
    return ontolist;
}

Rcpp::List sem1R::checkColDescription()
{
    if(printAndCheckPipeline(SET_DESCRIPTION_BICS) == PIPELINE_UNACCEPTED)
        return Rcpp::List();

    Rcpp::List ontolist(this->colOntology.size());
    int ionto = 0;
    for(boost::unordered_map<std::string, Ontology *>::iterator it = this->colOntology.begin(); it != this->colOntology.end(); ++it)
    {
        //refOntologies.push_back(it->second);

        Ontology *actOnto = it->second;
        std::vector<std::vector<std::string> > ontodescterms = actOnto->getDescriptionTerms();
        Rcpp::List listdesc(ontodescterms.size());
        for(int idesc = 0; idesc < ontodescterms.size(); ++idesc)
        {
            Rcpp::CharacterVector coldesc;
            for(int isubdesc = 0; isubdesc < ontodescterms[idesc].size(); ++isubdesc)
            {
                //if(actOnto->getOntologyParser()->getNodesBottomUP(&(ontodescterms[idesc][isubdesc])) != NULL)
                //Rcpp::Rcout << "#"  << ontodescterms[idesc][isubdesc] << "#- ";
                if(actOnto->getOntologyParser()->getNodesBottomUP(&(ontodescterms[idesc][isubdesc])) != NULL)
                {
                    //Rcpp::Rcout << "(" << ontodescterms[idesc][isubdesc] << " level: " << actOnto->getLevelOfSpecialization(actOnto->getOntologyParser()->getNodesBottomUP(&(ontodescterms[idesc][isubdesc]))) << ") ";
                    coldesc.push_back(ontodescterms[idesc][isubdesc] + std::string(", specificity level: ") + boost::lexical_cast<string>(actOnto->getLevelOfSpecialization(actOnto->getOntologyParser()->getNodesBottomUP(&(ontodescterms[idesc][isubdesc])))));
                }
                else
                {
                    //Rcpp::Rcout << "ontodesc: " << ontodescterms[idesc][isubdesc] << "not found!";
                    coldesc.push_back(ontodescterms[idesc][isubdesc] + std::string(" not found!"));
                }
            }
            listdesc[idesc] = coldesc;
            Rcpp::Rcout << std::endl;
        }
        ontolist[ionto] = listdesc;
        ++ionto;
    }
    return ontolist;
}

Rcpp::List sem1R::findDescription()
{
    if(printAndCheckPipeline(SET_DESCRIPTION_BICS) == PIPELINE_UNACCEPTED)
        return Rcpp::List();

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
    std::vector<Rcpp::CharacterVector> coveredNegRules(this->nrules);
    std::vector<Rcpp::DataFrame> igraphNodeVector;
    std::vector<Rcpp::DataFrame> igraphEdgeVector;
    
    Rcpp::CharacterVector dataRownames;
    Rcpp::CharacterVector dataColnames;

    if(!Rf_isNull(Rcpp::rownames(this->data)) || !Rf_isNull(Rcpp::colnames(this->data)))
    {
            dataRownames = Rcpp::rownames(this->data);
            dataColnames = Rcpp::colnames(this->data);
    }
    else
    {
        if(Rf_isNull(Rcpp::rownames(this->data)))
        {
            Rcpp::Rcerr << "Rownames of the input dataset are not set!" << std::endl;
        }
        if(Rf_isNull(Rcpp::colnames(this->data)))
        {
            Rcpp::Rcerr << "Rownames of the input dataset are not set!" << std::endl;
        }
        return Rcpp::List();
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
        bottomRules = bicExamples.initBottomFeaturesbitset(&enrichNodes, &(this->enrichPropositionTable3), 1);
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
        /*
        if(this->exhaustiveTest)
        {
                BeamTopDown *beamTDExh = new BeamTopDown(&bottomRules, &newclass, this->objective, this->verbose);
                if(this->verbose)
                    Rcpp::Rcout << std::endl << "################## EXHAUSTIVE VERSION ##################" << std::endl;
                beginExh = std::clock();
                newComplex actbestExh = beamTDExh->runExhaustive(this->filterTh, ruleDepth, this->signTH);
                endExh = std::clock();
        }
        */
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
        Rcpp::CharacterVector coveredNegList;
        if(this->ruleFormat == "col")
        {
            coveredList = beamTD->getOnlyCoveredColExamples(actbest, &newclass, &(armaDataTMP), dataColnames);
            coveredNegList = beamTD->getOnlyCoveredColNegExamples(actbest, &newclass, &(armaDataTMP), dataColnames);
        }
        else if(this->ruleFormat == "row")
        {
            coveredList = beamTD->getOnlyCoveredRowExamples(actbest, &newclass, &(armaDataTMP), dataRownames);
            coveredNegList = beamTD->getOnlyCoveredRowNegExamples(actbest, &newclass, &(armaDataTMP), dataRownames);
        }
        else if(this->ruleFormat == "both")
        {
            coveredList = beamTD->getOnlyCoveredRowColExamples(actbest, &newclass, &(armaDataTMP), dataRownames, dataColnames);
            coveredNegList = beamTD->getOnlyCoveredRowColNegExamples(actbest, &newclass, &(armaDataTMP), dataRownames, dataColnames);
        }
        else
            ;


        //igraph
        boost::unordered_map<std::string, boost::dynamic_bitset<> > ontologyCommonTerms;
        boost::unordered_map<std::string, boost::dynamic_bitset<> > ontologyRuleTerms;
        for(int iterm = 0; iterm < this->ruleset[irule].rules.size(); ++iterm)
        {           
            std::string ontoname = this->ruleset[irule].rules[iterm].noderef->onto_ref->getName();         
            if(ontologyCommonTerms.find(ontoname) == ontologyCommonTerms.end())
            {
                ontologyCommonTerms[ontoname] = boost::dynamic_bitset<>(this->ruleset[irule].rules[iterm].allGeneral.size());
                ontologyRuleTerms[ontoname] = boost::dynamic_bitset<>(this->ruleset[irule].rules[iterm].allGeneral.size());
            }
            //std::cout << " onto prev: " << ontologyCommonTerms[ontoname].count() << " new: " << this->ruleset[irule].rules[iterm].allGeneral.count() << std::endl;
            ontologyCommonTerms[ontoname] |= this->ruleset[irule].rules[iterm].allGeneral;
            ontologyCommonTerms[ontoname][this->ruleset[irule].rules[iterm].IDind] = 1;
            ontologyRuleTerms[ontoname][this->ruleset[irule].rules[iterm].IDind] = 1;
        }
        //iterate over commonterms
        Rcpp::CharacterVector termid;
        Rcpp::CharacterVector termname;
        Rcpp::CharacterVector edgefrom;
        Rcpp::CharacterVector edgeto;
        Rcpp::NumericVector termlevel;
        Rcpp::NumericVector isRuleterm;

        boost::unordered_map<std::string, boost::dynamic_bitset<> >::iterator itcommon;
        for(itcommon = ontologyCommonTerms.begin(); itcommon != ontologyCommonTerms.end();++itcommon)
        {
            //std::cout << "ONTO: " << itcommon->first << " generals: " << itcommon->second.count() <<  std::endl;
            size_t icombit = itcommon->second.find_first();
            while(icombit != boost::dynamic_bitset<>::npos)
            {
                //nodes
                termid.push_back(bottomRules[icombit].noderef->id);
                termname.push_back(bottomRules[icombit].noderef->name);
                termlevel.push_back(bottomRules[icombit].level);
                if(ontologyRuleTerms[itcommon->first][icombit])
                    isRuleterm.push_back(1);
                else
                    isRuleterm.push_back(0);
                //edges
                size_t iparent = bottomRules[icombit].allParents.find_first();
                while(iparent != boost::dynamic_bitset<>::npos)
                {
                    edgeto.push_back(bottomRules[iparent].noderef->id);
                    edgefrom.push_back(bottomRules[icombit].noderef->id);
                    //next
                    iparent = bottomRules[icombit].allParents.find_next(iparent);
                }
                //next
                icombit = itcommon->second.find_next(icombit);
            }
        }
        Rcpp::DataFrame nodes = Rcpp::DataFrame::create(Rcpp::Named("nodeID") = termid, Rcpp::_["nodeName"] = termname, Rcpp::_["nodeLevel"] = termlevel, Rcpp::_["isRuleTerm"] = isRuleterm);
        Rcpp::DataFrame edges = Rcpp::DataFrame::create(Rcpp::Named("edgeFrom") = edgefrom, Rcpp::_["edgeTo"] = edgeto);
        igraphEdgeVector.push_back(edges);
        igraphNodeVector.push_back(nodes);

        coveredRules[irule] = coveredList;
        coveredNegRules[irule] = coveredNegList;
        coveredString[irule] = coveredDesc;
					
        arma::mat new_armaData = beamTD->removeCoveredExamples(actbest, &(armaDataTMP), newclass);

        armaDataTMP = new_armaData;
        
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        //Rcpp::Rcout << "[TOPOLOGY VERSION] elapsed overall time: " << elapsed_secs;

        /*
        if(this->exhaustiveTest)
        {
                double elapsed_secsExh = double(endExh - beginExh) / CLOCKS_PER_SEC;
                //Rcpp::Rcout << std::endl << "[EXHAUSTIVE VERSION] elapsed overall time: " << elapsed_secsExh;
        }
        */
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
    Rcpp::List igraphNode(trueHypothesisiSize);
    Rcpp::List igraphEdge(trueHypothesisiSize);
    int ihypo = 0;
    for(int irule = 0; irule < trueHypothesisiSize; ++irule)
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
                Rcpp::_["nCoveredPOS"] = this->ruleset[irule].coverPos, Rcpp::_["nCoveredNEG"] = this->ruleset[irule].coverNeg,
                Rcpp::_["rules"] = getRuleID(&(this->ruleset[irule])), Rcpp::_["details"] = getRuleDetail(&(this->ruleset[irule])), Rcpp::_["coveredPOS"] = coveredRules[irule], Rcpp::_["coveredNEG"] = coveredNegRules[irule]);
                hypothesis[ihypo] = listRule;

                //igraph
                igraphNode[ihypo] = igraphNodeVector[ihypo];
                igraphEdge[ihypo] = igraphEdgeVector[ihypo];
                ++ihypo;
    }

    this->igraphEdge = igraphEdge;
    this->igraphNode = igraphNode;
    return hypothesis;
}
