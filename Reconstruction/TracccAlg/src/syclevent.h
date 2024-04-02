struct my_cell {
    uint64_t geometry_id = 0;
    uint64_t hit_id = 0;
    uint32_t channel0 = 0;
    uint32_t channel1 = 0;
    float timestamp = 0.;
    float value = 0.;
    uint64_t particle_id = 0;
};

class CEPCAlg{
public:
    virtual void initial(std::string* p_detect,
                         std::string* p_digit){}
    virtual void run(std::vector<my_cell> *cells, int cur_event){}
};

CEPCAlg* factory();