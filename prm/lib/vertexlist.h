#ifndef VERTEXLIST_H
#define VERTEXLIST_H

template<int ndof>
class vertexlist{
  struct node{
    //!q: array of structs
    float* q;
    int* env;
    node* next;
    node* parent;
  };



};

#endif // VERTEXLIST_H
