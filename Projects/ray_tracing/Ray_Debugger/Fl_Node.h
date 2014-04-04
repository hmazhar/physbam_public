#ifndef Fl_Node_H
#define Fl_Node_H

class Fl_Node {

  friend class Fl_Tree;

public:

  Fl_Node() {
    prev_ = 0;
    next_ = 0;
    sub_ = 0;
    vsub_ = 0;
    up_ = 0;
    opened_ = 0;
  }

protected:

  Fl_Node* prev_;
  Fl_Node* next_;
  Fl_Node* sub_;
  Fl_Node* vsub_;
  Fl_Node* up_;

  int opened_;

};

#endif
