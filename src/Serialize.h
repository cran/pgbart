#ifndef PGBART_SERIALIZE_H
#define PGBART_SERIALIZE_H

#include <cstdio>
#include <vector>

#include "Node.h"
#include "global.h"

class BartNode
{
public:
    int isbot;
    int var;
    double split;
    double node_mu;
    BartNode* LeftC;
    BartNode* RightC;

    bool Right(double *x){
        int i;
        if(x[this->var] > this->split)
            return true;
        else
            return false;
    }
};


void write_node(FILE* fp, Node* root){
    fprintf(fp, "%d ", root->Bot);
    if(root->Bot){
        fprintf(fp, "%f ", root->node_mu);
    }else{
        fprintf(fp, "%d %f ", root->rule.Var, RuleMat[root->rule.Var][root->rule.OrdRule]);
        write_node(fp, root->LeftC);
        write_node(fp, root->RightC);
    }
}


void add_itr(FILE* fp, std::vector<Node*> trees, int Ntrees){
    fprintf(fp, "\n");
    for(int i = 1; i <= Ntrees; i++){
        write_node(fp, trees[i]);
        fprintf(fp, "\n");
    }
}

void read_node(FILE* fp, BartNode* root){
    fscanf(fp, "%d", &(root->isbot));
    if(root->isbot){
        fscanf(fp, "%f", &(root->node_mu));
    }else{
        fscanf(fp, "%d %f", &(root->var), &(root->split));
        BartNode* left = new BartNode;
        BartNode* right = new BartNode;
        root->LeftC = left;
        root->RightC = right;
        read_node(fp, left);
        read_node(fp, right);
    }
}


void read_itr(FILE* fp, std::vector<BartNode*> trees, int Ntrees){
    int isbot, var;
    double split, node_mu;
    for(int i = 1; i <= Ntrees; i++){
        read_node(fp, trees[i]);
    }
}

void clear_node(BartNode* root){
    if(!(root->isbot)){
        clear_node(root->LeftC);
        clear_node(root->RightC);
    }
    delete root;
}

void clear_itr(std::vector<BartNode*> trees, int Ntrees){
    for(int i = 1; i <= Ntrees; i++)
        clear_node(trees[i]);
}

double predict(BartNode* node, double* x){
    if(node->isbot)
        return node->node_mu;
    else{
        if(node->Right(x))
            return predict(node->RightC, x);
        else
            return predict(node->LeftC, x);
    }
}

#endif
