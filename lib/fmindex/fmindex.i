/* File: fmindex.i */


%include stdint.i

%include std_vector.i
namespace std {
    %template(IntVector) vector<uint32_t>;
}

%include typemaps.i
%apply uint32_t *OUTPUT { uint32_t *matches };
%apply uint32_t *OUTPUT { uint32_t *max_suffix_length };

%module fmindex
%{
#define SWIG_FILE_WITH_INIT
class FM {
public:
    FM(uint8_t* T,uint32_t n,uint32_t samplerate);
    void build(uint8_t* T,uint32_t n,uint32_t samplerate);
    static FM* load(char* filename);
    int32_t save(char* filename);
    uint8_t* remap0(uint8_t* T,uint32_t n);
    uint32_t count(char* pattern);
    uint32_t count_max_suffix(char* pattern, uint32_t* max_suffix_length);	
    std::vector<uint32_t> locate(char* pattern,uint32_t* matches);
    uint32_t getSize();
    char* extract(uint32_t start,uint32_t stop);
    uint8_t* reconstructText(uint32_t* n);
    float getSizeN();
    virtual ~FM();
};

%}

class FM {
public:
    FM(uint8_t* T,uint32_t n,uint32_t samplerate);
    void build(uint8_t* T,uint32_t n,uint32_t samplerate);
    static FM* load(char* filename);
    int32_t save(char* filename);
    uint8_t* remap0(uint8_t* T,uint32_t n);
    uint32_t count(char* pattern);
    uint32_t count_max_suffix(char* pattern, uint32_t* max_suffix_length);
    std::vector<uint32_t>  locate(char* pattern,uint32_t* matches);
    uint32_t getSize();
    char* extract(uint32_t start,uint32_t stop);
    uint8_t* reconstructText(uint32_t* n);
    float getSizeN();
    virtual ~FM();
};
