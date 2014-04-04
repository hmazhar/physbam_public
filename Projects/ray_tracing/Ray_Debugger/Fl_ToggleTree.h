#ifndef Fl_ToggleTree_H
#define Fl_ToggleTree_H

#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include "Fl_ToggleNode.h"
#include "Fl_Tree.h"
#include <FL/Fl_Pixmap.H>

enum Fl_ToggleState {
  FL_TOGGLE_NONE = 0,
  FL_TOGGLE_SELECT = 1,
  FL_TOGGLE_RESELECT = 2,
  FL_TOGGLE_SELECT_MASK = 3,
  FL_TOGGLE_OPENED = 4,
  FL_TOGGLE_CLOSED = 8,
  FL_TOGGLE_HIT = 16
};

class Fl_Input;

class Fl_ToggleTree : public Fl_Tree {

public:

  Fl_ToggleTree(int x, int y, int w, int h,const char *label=0);
  virtual ~Fl_ToggleTree();

  int handle(int event) PHYSBAM_OVERRIDE;

  void resize(int x,int y,int w,int h)
  {
//    printf("%d %d %d %d\n",x,y,w,h);
      Fl_Tree::resize(x,y,w,h);
  }

  Fl_ToggleNode* current(void) {
    return (Fl_ToggleNode *)t_current_;
  }
  Fl_ToggleState state(void) {
    return state_;
  }

  void selection_label_color(Fl_Color c) {
    selection_label_color_ = c;
  }
  Fl_Color selection_label_color(void) {
    return selection_label_color_;
  }
  void alternate_color(Fl_Color c) {
    alternate_color_ = c;
  }
  Fl_Color alternate_color(void) {
    return alternate_color_;
  }

  void trim_color(Fl_Color c) {
    trim_color_ = c;
  }
  Fl_Color trim_color(void) {
    return trim_color_;
  }

  void indent_toggles(int b) {
    indent_toggles_ = b;
  };
  int indent_toggles(void) {
    return indent_toggles_;
  };

  void draw_lines(int b) {
    draw_lines_ = b;
  };
  int draw_lines(void) {
    return draw_lines_;
  };

  void open(Fl_ToggleNode* node);
  void close(Fl_ToggleNode* node);

  void label_offset(int l) {
    label_offset_ = l;
  }
  void pixmap_offset(int l) {
    pixmap_offset_ = l;
  }
  void opened_pixmap(Fl_Pixmap *);
  void closed_pixmap(Fl_Pixmap *a);

  Fl_Pixmap * opened_pixmap() {
    return opened_pixmap_;
  }
  Fl_Pixmap * closed_pixmap() {
    return closed_pixmap_;
  }

  const int* column_widths() const {return column_widths_; }
  void column_widths(const int* l) {
    column_widths_ = l;
  }

  Fl_Font textfont() const {return (Fl_Font)textfont_; }
  void textfont(uchar s) {
    textfont_ = s;
  }
  uchar textsize() const {return textsize_; }
  void textsize(uchar s) {
    textsize_ = s;
  }
  Fl_Color textcolor() const {return (Fl_Color)textcolor_; }
  void textcolor(uchar n) {
    textcolor_ = n;
  }

  char column_char() const {return column_char_; }
  void column_char(char c) {
    column_char_ = c;
  }

  Fl_ToggleNode* selection(void);
  Fl_ToggleNode* selection(int i);
  int selection_count(void);

  Fl_ToggleNode* selected(void) {
    return (Fl_ToggleNode *)current_;
  }
  void select_range(Fl_ToggleNode* start, Fl_ToggleNode* end, int add = 0);
  void select(Fl_ToggleNode* s) {
    select_range(s,s,0);
  }
  void unselect(void) {
    select_range(0, 0, 0);
  }

  /* Remove an item from the tree but also select the previous line
     if possible. */

  int remove (Fl_ToggleNode * a, bool dodelete=1) {
    Fl_ToggleNode * sel = 0;
    if (a->selected_) {
      if (a->up_)
    sel = (Fl_ToggleNode *)a->up_;
      else if (a->prev_)
    sel = (Fl_ToggleNode *)a->prev_;
    }
    Fl_Tree::remove((Fl_Node *) a,dodelete);

    if (!sel) sel = (Fl_ToggleNode *)first_;
    if (sel) sel->selected_ = 1;
    current_ = sel;
    redraw();
    return 1;
  }
  int remove (void * a);
  int remove (char * a);

  void traverse_start(Fl_Node * a) {
    Fl_Tree::traverse_start(a);
  }
  Fl_ToggleNode * traverse_start() {
    return (Fl_ToggleNode *) Fl_Tree::traverse_start();
  }
  Fl_ToggleNode * traverse_forward(int visible, int &depth) {
    return (Fl_ToggleNode *) Fl_Tree::traverse_forward(visible, depth);
  }
  Fl_ToggleNode * traverse_forward() {
    return (Fl_ToggleNode *) Fl_Tree::traverse_forward();
  }
  Fl_ToggleNode * traverse_backward() {
    return (Fl_ToggleNode *) Fl_Tree::traverse_backward();
  }

  Fl_ToggleNode * find (void * a);
  Fl_ToggleNode * find (char * a);

  Fl_ToggleNode * add_sub(char* label = 0, int can_open = 1,
                          Fl_Pixmap* pixmap = 0, void * d = 0) {
    Fl_ToggleNode * node;
    Fl_Tree::add_sub(node = new Fl_ToggleNode(label, can_open, pixmap, d));
    return node;
  }

  Fl_ToggleNode * add_next(char* label = 0, int can_open = 1,
                           Fl_Pixmap* pixmap = 0, void * d = 0) {
    Fl_ToggleNode * node;
    Fl_Tree::add_next(node = new Fl_ToggleNode(label, can_open, pixmap, d));
    return node;
  }

  void add_sub(Fl_ToggleNode * a) {
    Fl_Tree::add_sub(a);
  }

  void add_next(Fl_ToggleNode * a) {
    Fl_Tree::add_next(a);
  }

  static int sort_by_label(Fl_Node* a, Fl_Node* b);

  void edit_on_reselect(int b) {
    edit_on_reselect_ = b;
  }

  void edit_callback(Fl_Callback* c, void* p);
  void edit_callback(Fl_Callback* c);
  void edit_callback(Fl_Callback0*c);
  void edit_callback(Fl_Callback1*c, long p = 0);

  static void edit_default_callback(Fl_Input* input, void* ptr);
  void end_edit(void);

protected:
  virtual void edit(Fl_ToggleNode* t, int cx, int cy);

  void draw_label(char* str, int indent, int x, int y, int w, int h);
  virtual void draw_node(int depth, int cy, Fl_Node* node);
  int label_offset_;
  int pixmap_offset_;
  int indent_toggles_;
  int edit_on_reselect_;
  int draw_lines_;
  Fl_Color alternate_color_;
  Fl_Color selection_label_color_;
  Fl_Color trim_color_;
  const int* column_widths_;
  char column_char_;

  uchar textfont_, textsize_, textcolor_;

  Fl_ToggleState state_;
  Fl_Pixmap* opened_pixmap_;
  Fl_Pixmap* closed_pixmap_;
  static Fl_Pixmap* s_opened_pixmap_;
  static Fl_Pixmap* s_closed_pixmap_;

  int selection_i_;
  int selection_count_;
  Fl_ToggleNode* selection_current_;

  Fl_Input *edit_input_;
};

#endif
