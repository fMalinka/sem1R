#include "error.h"

int pipelineCODE[5] = {0,0,0,0,0};
int printAndCheckPipeline(int actState)
{
    if(actState == SET_DATASET)
    {
        for(int istat = 0; istat < 5; ++istat)
        {
            pipelineCODE[istat] = 0;
        }
    }
    else if(actState == SET_COL_ONTOLOGY)
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
    else if(pipelineCODE[SET_DESCRIPTION_BICS] == 0)
    {
        if(pipelineCODE[SET_DATASET] == 0)
        {
            Rcpp::Rcerr << "Dataset must be loaded first! Use setDataset function." << std::endl;
            return PIPELINE_UNACCEPTED;
        }

        if(pipelineCODE[SET_COL_ONTOLOGY] == 0 && pipelineCODE[SET_ROW_ONTOLOGY] == 0)
        {
            Rcpp::Rcerr << "At least column or row ontology must be created first! Use createCOLOntology or createROWntology function!" << std::endl;
            return PIPELINE_UNACCEPTED;
        }
    }

    pipelineCODE[actState] = 1;
    return PIPELINE_OK;
}

