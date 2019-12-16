#include "Evaluate.h"

/*
 * The code has been written for minimization of objectives.
 * For maximalization use negative of the function value as the obj. value.
 */

Evaluate::Evaluate()
{
    ;
}

//Return size of bicluster given by individual
int Evaluate::evaluateSize(arma::vec *gene, arma::vec *sample)
{
/*    int sumrow = 0;
    int sumcol = 0;
    //iterate over rows/genes
    for(int bitrow = 0; bitrow < nrow; ++bitrow)
    {
        if(gene[bitrow] == 1)
            ++sumrow;
    }
    //iterate over columns/samples
    for(int bitcol = 0; bitcol < ncol; ++bitcol)
    {
        if(sample[bitcol] == 1)
            ++sumcol;
    }
   */
    //return -(sumrow*sumcol);
    return -(arma::accu(*gene) * arma::accu(*sample));
}

//Return Number of correct predictions
int Evaluate::evaluateCorrectPrediction(arma::vec *gene, arma::vec *sample, arma::mat *data)
{
/*    arma::vec rowBitset(nrow, arma::fill::zeros);
    arma::vec colBitset(ncol, arma::fill::zeros);
    bool anyRow = false;
    bool anyCol = false;
    //iterate over rows/genes
    for(int bitrow = 0; bitrow < nrow; ++bitrow)
    {
        //for more speed use [] instead of () -> bound check will be disable
        if(gene[bitrow] == 1)
        {
            rowBitset(bitrow) = 1;
            anyRow = true;
        }
    }
    //iterate over columns/samples
    for(int bitcol = 0; bitcol < ncol; ++bitcol)
    {
        //for more speed use [] instead of () -> bound check will be disable
        if(sample[bitcol] == 1)
        {
            colBitset(bitcol) = 1;
            anyCol = true;
        }
    }

    //invalid bicluster form
    if(!anyCol || !anyRow)
    {
        //std::cout << "PRAZDNO!!" << std::endl;
        return 1;
    }
*/
    if(arma::accu(*gene) == 0 ||  arma::accu(*sample) == 0)
        return 1;
    arma::mat binTemplate = (*gene) * arma::trans(*sample);  //define matrix representing of selected elements
    arma::mat correctPredictData = binTemplate % *data; //% indicates element-wise multiplication

    return -(arma::accu(correctPredictData));
}

//Return Number of correct predictions - for SIGNIFICANT TEST
int Evaluate::evaluateCorrectPrediction(Rcpp::NumericVector *tmpGene, Rcpp::NumericVector *tmpSample, arma::mat *data)
{
    arma::vec rowBitset = (arma::vec)*tmpGene;
    arma::vec colBitset = (arma::vec)*tmpSample;

    arma::mat binTemplate = rowBitset * arma::trans(colBitset);  //define matrix representing of selected elements
    //arma::mat originalData = as<arma::mat>(*data);
    arma::mat correctPredictData = binTemplate % *data; //% indicates element-wise multiplication

    return -(arma::accu(correctPredictData));
}

//Return Number of different predictions
int Evaluate::evaluateXOR(int *gene, int *sample, int nrow, int ncol, arma::mat *data)
{
    arma::vec rowBitset(nrow, arma::fill::zeros);
    arma::vec colBitset(ncol, arma::fill::zeros);
    //iterate over rows/genes
    for(int bitrow = 0; bitrow < nrow; ++bitrow)
    {
        //for more speed use [] instead of () -> bound check will be disable
        if(gene[bitrow] == 1)
            rowBitset(bitrow) = 1;
    }
    //iterate over columns/samples
    for(int bitcol = 0; bitcol < ncol; ++bitcol)
    {
        //for more speed use [] instead of () -> bound check will be disable
        if(sample[bitcol] == 1)
            colBitset(bitcol) = 1;
    }

    arma::mat binTemplate = rowBitset * arma::trans(colBitset);  //define matrix representing of selected elements
    //arma::mat originalData = as<arma::mat>(*data);

    arma::uvec indeces = find(binTemplate); //find non zero --- indexes
    arma::vec elementsTempl = (binTemplate)(indeces);

    arma::vec origVector = (arma::vectorise(*data));
    arma::vec elementsOrig = origVector(indeces);  //vector with selected elements

    return arma::accu(arma::abs(elementsOrig-elementsTempl));
}

//Return
double Evaluate::evaluateF1score(arma::vec *gene, arma::vec *sample, arma::mat *data)
{
/*    arma::vec rowBitset(nrow, arma::fill::zeros);
    arma::vec colBitset(ncol, arma::fill::zeros);
    //iterate over rows/genes
    for(int bitrow = 0; bitrow < nrow; ++bitrow)
    {
        //for more speed use [] instead of () -> bound check will be disable
        if(gene[bitrow] == 1)
            rowBitset(bitrow) = 1;
    }
    //iterate over columns/samples
    for(int bitcol = 0; bitcol < ncol; ++bitcol)
    {
        //for more speed use [] instead of () -> bound check will be disable
        if(sample[bitcol] == 1)
            colBitset(bitcol) = 1;
    }
*/
    arma::mat binTemplate = (*gene) * arma::trans(*sample);  //define matrix representing of selected elements
    //arma::mat originalData = as<arma::mat>(*data);

    double TP = arma::accu(binTemplate % *data);
    arma::vec origVector = (arma::vectorise(*data));
    arma::vec templVector = (arma::vectorise(binTemplate));
    arma::vec xr = origVector - templVector;
    //xr = xr % templVector; //elimate elements outsite of bicluster
    //double FN = arma::size(arma::find(xr > 0));
    //double FP = arma::size(arma::find(xr < 0));
    arma::uvec vecFN = arma::find(xr > 0);
    arma::uvec vecFP = arma::find(xr < 0);
    double FN = vecFN.n_elem;
    double FP = vecFP.n_elem;
    
    double precision = TP/(TP+FP);
    double recall = TP/(TP+FN);

    //std::cout << "TP: " << TP << "FN: " << FN << "FP: " << FP << " prec: " << precision << " re: " << recall << "f1: " << -2*((precision*recall)/(precision + recall)) << std::endl;

    return -2*((precision*recall)/(precision + recall));
}

//Return Number of different predictions - for SIGNIFICANT TEST
int Evaluate::evaluateXOR(Rcpp::NumericVector *tmpGene, Rcpp::NumericVector *tmpSample, arma::mat *data)
{
    arma::vec rowBitset = (arma::vec)*tmpGene;
    arma::vec colBitset = (arma::vec)*tmpSample;

    arma::mat binTemplate = rowBitset * arma::trans(colBitset);  //define matrix representing of selected elements
    //arma::mat originalData = as<arma::mat>(*data);

    arma::uvec indeces = find(binTemplate); //find non zero --- indexes
    arma::vec elementsTempl = (binTemplate)(indeces);

    arma::vec origVector = (arma::vectorise(*data));
    arma::vec elementsOrig = origVector(indeces);  //vector with selected elements

    return arma::accu(arma::abs(elementsOrig-elementsTempl));
}


//Return Semantic of bicluster given by individual
double Evaluate::evaluateSemantic(arma::vec *gene, int nStart, int nStop, boost::unordered_map<std::string, Ontology *> *ontology)
{
    std::vector<double> score;   //final score

    for(boost::unordered_map<std::string, Ontology *>::iterator it = ontology->begin(); it != ontology->end(); ++it)
    {
        Ontology *onto = it->second;

        //iterate over rows/genes
        int ibit = 0;
        std::vector<int> bit;
        int i = 0;
        for(int bitrow = nStart; bitrow < nStop; ++bitrow)
        {
            if((*gene)[bitrow] == 1)
                bit.push_back(i);
            i++;
        }

        double actscore = 0;
        //jak se zachovat pri vyberu 1 sloupec/radku?
        for(int r = 0; r < bit.size(); ++r)
        {
            for(int c = r+1; c < bit.size(); ++c)
            {
                //spatne indexujes, mel bys indexovat pozici (bitrow)
                actscore += onto->getSemanticDistance(bit[r],bit[c]);
                ibit++;
            }
        }
        if(ibit == 0)
        {
            //std::cout << "NULA v1" << std::endl;
            score.push_back(0);
        }
        else
            score.push_back(actscore/(double)ibit); //average score for the ontology
    }

    double finalScore = 0;
    for(int i = 0; i < score.size(); ++i)
    {
        finalScore += score[i];
    }
    //std::cout << "score: " << -(finalScore/(double)score.size()) <<std::endl;
    if(score.size() == 0)
    {
        //std::cout << "NULA v2" << std::endl;
        return 0;
    }

    return -(finalScore/(double)score.size());   //total score divided by number of ontologies
}

//Return Semantic of bicluster given by individual - for SIGNIFICANT TEST
double Evaluate::evaluateSemantic(Rcpp::NumericVector *bitset, boost::unordered_map<std::string, Ontology *> *ontology)
{
    std::vector<double> score;   //final score

    for(boost::unordered_map<std::string, Ontology *>::iterator it = ontology->begin(); it != ontology->end(); ++it)
    {
        Ontology *onto = it->second;

        //iterate over rows/genes
        int ibit = 0;
        std::vector<int> bit;
        for(int bitrow = 0; bitrow < bitset->length(); ++bitrow)
        {
            if((*bitset)[bitrow] == 1)
                bit.push_back(bitrow);
        }

        double actscore = 0;
        for(int r = 0; r < bit.size(); ++r)
        {
            for(int c = r+1; c < bit.size(); ++c)
            {
                actscore += onto->getSemanticDistance(r,c);
                ibit++;
            }
        }
        if(ibit == 0)
            score.push_back(0);
        else
            score.push_back(actscore/(double)ibit); //average score for the ontology
    }

    double finalScore = 0;
    for(int i = 0; i < score.size(); ++i)
    {
        finalScore += score[i];
    }
    //std::cout << "score: " << -(finalScore/(double)score.size()) <<std::endl;
    if(score.size() == 0)
        return 0;

    return -(finalScore/(double)score.size());   //total score divided by number of ontologies
}


void Evaluate::updataDataMatrix(arma::vec *gene, arma::vec *sample, arma::mat *data)
{
 /*   arma::vec rowBitset(nrow, arma::fill::zeros);
    arma::vec colBitset(ncol, arma::fill::zeros);
    //iterate over rows/genes
    for(int bitrow = 0; bitrow < nrow; ++bitrow)
    {
        //for more speed use [] instead of () -> bound check will be disable
        if(gene[bitrow] == 1)
            rowBitset(bitrow) = 1;
    }
    //iterate over columns/samples
    for(int bitcol = 0; bitcol < ncol; ++bitcol)
    {
        //for more speed use [] instead of () -> bound check will be disable
        if(sample[bitcol] == 1)
            colBitset(bitcol) = 1;
    }
*/
    arma::mat binTemplate = (*gene) * arma::trans(*sample);  //define matrix representing of selected elements

    arma::mat commonElements = *data % binTemplate;
    arma::mat matOnes(data->n_rows, data->n_cols, arma::fill::ones);
    arma::mat NEGcommonElements = arma::abs(commonElements - matOnes);
    *data = *data % NEGcommonElements;
}
