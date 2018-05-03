//=======================================================================================
//  This is enum class definition for unit test.
//=======================================================================================

enum class TGT_STATE {
    A,
    B,
    C,
    D,
};
namespace ENUM {

    static const std::unordered_map<std::string, TGT_STATE> table_str_TGT_STATE{
        {"A", TGT_STATE::A},
        {"B", TGT_STATE::B},
        {"C", TGT_STATE::C},
        {"D", TGT_STATE::D},
    };

    static const std::unordered_map<TGT_STATE, std::string> table_TGT_STATE_str{
        {TGT_STATE::A, "A"},
        {TGT_STATE::B, "B"},
        {TGT_STATE::C, "C"},
        {TGT_STATE::D, "D"},
    };

    TGT_STATE which_TGT_STATE(const std::string &str){
        if(table_str_TGT_STATE.find(str) != table_str_TGT_STATE.end()){
            return table_str_TGT_STATE.at(str);
        } else {
            std::cerr << "  TGT_STATE: input = " << str << std::endl;
            throw std::out_of_range("undefined enum value in TGT_STATE.");
        }
    }

    std::string what(const TGT_STATE &e){
        if(table_TGT_STATE_str.find(e) != table_TGT_STATE_str.end()){
            return table_TGT_STATE_str.at(e);
        } else {
            using type_base = typename std::underlying_type<TGT_STATE>::type;
            std::cerr << "  TGT_STATE: input = " << static_cast<type_base>(e) << std::endl;
            throw std::out_of_range("undefined enum value in TGT_STATE.");
        }
    }
}
std::ostream& operator << (std::ostream &os, const TGT_STATE &v){
    os << ENUM::what(v);
    return os;
}
namespace std {
    inline string to_string(const TGT_STATE &e){ return ENUM::what(e); }
}
