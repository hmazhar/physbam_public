#ifndef __FLTTREE_RESIZER_H
#define __FLTTREE_RESIZER_H
class Fl_ToggleTreeResizer:public Fl_Widget
{
    Fl_Scroll* scroll_;
    Fl_ToggleTree* tt_;

    int xo;
    int minw;
public:
    Fl_ToggleTreeResizer(Fl_Scroll* scroll,Fl_ToggleTree* tt)
    : Fl_Widget(scroll->x(),scroll->y(),scroll->w(),scroll->h()) 
    {
        scroll_=scroll;
        tt_=tt;
        
        minw=tt_->w();
        xo=tt_->x()-scroll_->x();
    }
    void draw(){}
    void resize(int X,int Y,int W,int H) 
    {
        W=W-scroll_->scrollbar.w()-xo-1;
        if (W<minw) W=minw;
        tt_->size(W,tt_->h());
    }
};
#endif
