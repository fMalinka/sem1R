#include "error.h"

int pipelineCODE[5] = {0,0,0,0,0};
int printAndCheckPipeline(int actState)
{
    if(actState == SET_COL_ONTOLOGY)
    {
        if(pipelineCODE[SET_DATASET] == 0)
        {
            Rcpp::Rcerr << "Dataset must be loaded first! Use setDataset function." << std::endl;
            return PIPELINE_UNACCEPTED;
        }
    }
    else if(actState == SET_ROW_ONTOLOGY)
    {
        if(pipelineCODE[SET_DATASET] == 0)
        {
            Rcpp::Rcerr << "Dataset must be loaded first! Use setDataset function." << std::endl;
            return PIPELINE_UNACCEPTED;
        }
    }
    else if(actState == SET_FIND_BICS)
    {
        if(pipelineCODE[SET_DATASET] == 0)
        {
            Rcpp::Rcerr << "Dataset must be loaded first! Use setDataset function." << std::endl;
            return PIPELINE_UNACCEPTED;
        }

        if(pipelineCODE[SET_COL_ONTOLOGY] == 0)
        {
            Rcpp::Rcerr << "Column ontology must be created first! Use createCOLOntology function!" << std::endl;
            return PIPELINE_UNACCEPTED;
        }

        if(pipelineCODE[SET_ROW_ONTOLOGY] == 0)
        {
            Rcpp::Rcerr << "Row ontology must be created first! Use createROWOntology function!" << std::endl;
            return PIPELINE_UNACCEPTED;
        }
    }
    else if(pipelineCODE[SET_DESCRIPTION_BICS] == 0)
    {
        if(pipelineCODE[SET_DATASET] == 0)
        {
            Rcpp::Rcerr << "Dataset must be loaded first! Use setDataset function." << std::endl;
            return PIPELINE_UNACCEPTED;
        }

        if(pipelineCODE[SET_COL_ONTOLOGY] == 0)
        {
            Rcpp::Rcerr << "Column ontology must be created first! Use createCOLOntology function!" << std::endl;
            return PIPELINE_UNACCEPTED;
        }

        if(pipelineCODE[SET_ROW_ONTOLOGY] == 0)
        {
            Rcpp::Rcerr << "Row ontology must be created first! Use createROWOntology function!" << std::endl;
            return PIPELINE_UNACCEPTED;
        }

       if(pipelineCODE[SET_FIND_BICS] == 0)
       {
           Rcpp::Rcerr << "Pareto set must be computed first! Use findPareto function!" << std::endl;
           return PIPELINE_UNACCEPTED;
       }
    }

    pipelineCODE[actState] = 1;
    return PIPELINE_OK;
}

