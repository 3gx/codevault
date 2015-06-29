#pragma once
#include <sstream>
#include <tuple>
#include <functional>
#include <iostream>

namespace parse_arguments
{
  // -----------------------------
  //  Converts from string to type
  //    Should be specialized for user-defined types
  // -----------------------------
  template<typename T> 
    T string_to_type(std::string s) { return static_cast<T>(std::stod(s)); }
  template<> std::string string_to_type<std::string>(std::string s) { return s; };
  template<> bool string_to_type<bool>(std::string s)  { return true; };

  // ----------------
  //  Pretty printer
  //    Should be specialized for user-defined types
  // ----------------
  template<typename T>  
    std::string pretty_print(T t) { return std::to_string(t); }
  template<> std::string pretty_print<bool>(bool t) { return t ? "true" : "false"; }
  template<> std::string pretty_print<std::string>(std::string t) { return "'"+t+"'"; } 


  /********************************
   *  -= parser specific code =-  *
   ********************************/

  // ----------------
  //  Parameter name
  // ----------------
  class ParamName
  {
    private:
      std::string const description_;
      std::string const short_name_;
      std::string const long_name_;
      std::string arg_string_;
    public:
      ParamName(std::string desc, std::string nshort, std::string nlong) :
        description_(std::move(desc)),
        short_name_(std::move(nshort)),
        long_name_(std::move(nlong)) {}

      enum NameType {DESCRIPTION, SHORT_NAME, LONG_NAME, ARG_STRING};
      std::string get(NameType name) const 
      {
        switch (name)
        {
          case DESCRIPTION: return description_;
          case SHORT_NAME:  return short_name_;
          case LONG_NAME:   return long_name_;
          case ARG_STRING:  return arg_string_;
        }
        return std::string{};
      }
      std::string& getArgString()
      {
        return arg_string_;
      }
  };

  // --------------
  //  Parser class
  // --------------
  template<typename... Ts>
  class Parser
  {
    private:
      using ParamPack = std::tuple<Ts...>;
      static constexpr auto ParamSize = sizeof...(Ts);

      int argc_;
      char ** argv_;
      ParamPack param_pack_;

      // -----------------------------
      //  Actual argument parser code
      // -----------------------------
      template<typename T>
        static bool getParam(std::string arg, std::string what, T& value, std::string& parm)
        {
          if (arg.find(what) == std::string::npos)
            return false;
          value = string_to_type<T>(arg.substr(what.size())); 
          parm = what;
          parm.erase(what.size()-1);
          return true;
        }

      template<size_t I>
        auto parseImpl(std::string arg, bool flag) -> 
        typename std::enable_if<(I == ParamSize),bool>::type
        {
          return flag;
        }
      template<size_t I=0>
        auto parseImpl(std::string arg, bool flag = false) -> 
        typename std::enable_if<(I < ParamSize),bool>::type
        {
          using namespace std;
          auto& t = get<I>(param_pack_);
          auto const & shortParam = t.get(ParamName::SHORT_NAME);
          auto const & longParam = t.get(ParamName::LONG_NAME);

          using value_type = typename tuple_element<I,ParamPack>::type::value_type; 
          string ss_eq="=", arg1=arg;
          if (is_same<value_type,bool>::value)
          {
            arg1 += "\\";
            ss_eq = "\\";
          }
          if (!shortParam.empty())
            flag |= getParam(arg1,
                string("-")+shortParam+ss_eq,
                t.value(), t.getArgString());
          if (!longParam.empty())
            flag |= getParam(arg1, 
                string("--")+longParam+ss_eq,
                t.value(), t.getArgString());
          return parseImpl<I+1>(arg, flag);
        }
  
      // ------------------------------
      //  Find max argument name width
      // ------------------------------
      template<size_t I>
        auto maxWidthNameImpl(ParamName::NameType name, size_t count) ->
        typename std::enable_if<(I == ParamSize),size_t>::type
        {
          return count;
        }
      template<size_t I=0>
        auto maxWidthNameImpl(ParamName::NameType name, size_t count = 0) ->
        typename std::enable_if<(I < ParamSize),size_t>::type
        {
          count = std::max(count, std::get<I>(param_pack_).get(name).size());
          return maxWidthNameImpl<I+1>(name, count);
        }
      size_t maxWidthName(ParamName::NameType name)
      {
        return maxWidthNameImpl(name);
      }

      // ---------------
      //  Usage printer
      // ---------------
      template<size_t I, typename SS>
        auto printUsageImpl(SS&& ss, size_t maxWidthLong, size_t maxWidthShort) ->
        typename std::enable_if<(I == ParamSize)>::type
        {}
      template<size_t I=0, typename SS>
        auto printUsageImpl(SS&& ss, size_t maxWidthLong, size_t maxWidthShort) ->
        typename std::enable_if<(I < ParamSize)>::type
        {
          using namespace std;
          auto const & t = get<I>(param_pack_);
          using value_type = typename tuple_element<I,ParamPack>::type::value_type; 

          string ss_value;
          if (!is_same<value_type,bool>::value)
            ss_value = "=";
          else
            ss_value = " ";

          ss << "  ";
          if (!t.get(ParamName::SHORT_NAME).empty())
            ss << "-" << t.get(ParamName::SHORT_NAME) << ss_value;
          else
            ss << "  ";
          {
            auto width = t.get(ParamName::SHORT_NAME).size();
            while(width++ < maxWidthShort)
              ss << " ";
          }
          if (!t.get(ParamName::LONG_NAME).empty())
            ss << "  --" << t.get(ParamName::LONG_NAME) << ss_value;
          else 
            ss << "     ";

          {
            auto width = t.get(ParamName::LONG_NAME).size();
            while(width++ < maxWidthLong+4)
              ss << " ";
          }

          ss << t.get(ParamName::DESCRIPTION)
            << " [";
          ss << pretty_print(t.default_value());
          ss << "]\n";
          return printUsageImpl<I+1>(ss, maxWidthLong,maxWidthShort);
        }
      // -------------------
      //  Parameter printer
      // -------------------
      template<size_t I, typename SS>
        auto printParamsImpl(SS&& ss, size_t wD, size_t wA) ->
        typename std::enable_if<(I == ParamSize)>::type
        {}
      template<size_t I=0, typename SS>
        auto printParamsImpl(SS&& ss, size_t wD, size_t wA) ->
        typename std::enable_if<(I < ParamSize)>::type
        {
          using namespace std;
          auto const & t = get<I>(param_pack_);

          ss << "\t";
          ss << t.get(ParamName::DESCRIPTION) << ":";
          auto width = t.get(ParamName::DESCRIPTION).size();
          while(width++ < wD+2)
            ss << " ";
          
          using value_type = typename tuple_element<I,ParamPack>::type::value_type; 
          string ss_value = "", ss_space=" ";

          if (!is_same<value_type,bool>::value && !t.get(ParamName::ARG_STRING).empty())
          {
            ss_value = "=";
            ss_space = "";
          }

          {
            auto width = t.get(ParamName::ARG_STRING).size();
            while(width++ < wA)
              ss << " ";
          }
          ss << ss_space << t.get(ParamName::ARG_STRING) << ss_value;
          ss << "  ";

          ss << pretty_print(t.value());
          if (t.default_value() == t.value())
            ss << " [=default]";

          ss << endl;

          printParamsImpl<I+1>(ss, wD, wA);
        }

    public:

      Parser(int argc, char** argv, Ts&&... ts) :
        argc_{argc}, argv_{argv},
        param_pack_{std::make_tuple(std::forward<Ts>(ts)...)}
      {}

      bool parse()
      {
        using namespace std;
        for (int i = 1; i < argc_; ++i)
          if (!parseImpl(argv_[i]))
            return false;
        return true;
      };
  
      std::string usage()
      {
        using namespace std;
        stringstream ss;
        ss << endl;
        ss << "Usage: " ;
        for (int i = 0; i < argc_; ++i)
          ss << argv_[i] << " ";
        ss << endl;
        ss << endl;
        ss << "Available arguments: \n";
        auto const maxWidthLongName  = maxWidthName(ParamName::LONG_NAME);
        auto const maxWidthShortName = maxWidthName(ParamName::SHORT_NAME);
        printUsageImpl(ss,maxWidthLongName,maxWidthShortName);
        ss << endl;
        return ss.str();
      }
 
      template<typename F>
        std::string parse_all(F f)
        {
          if (!parse())
          {
            f(usage());
            return std::string{};
          }
          return params();
        }

      std::string parse_all()
      {
        return parse_all([](std::string s) { std::cerr << s; exit(-1); });
      }


      std::string params()
      {
        using namespace std;
        auto const maxWidthArgString   = maxWidthName(ParamName::ARG_STRING);
        auto const maxWidthDescription = maxWidthName(ParamName::DESCRIPTION);
        stringstream ss;
        ss << endl;
        ss << "Exec line: ";
        for (int i = 0; i < argc_; ++i)
          ss << argv_[i] << " ";
        ss << endl << endl;
        ss << "Parameters: \n";
        printParamsImpl(ss, maxWidthDescription, maxWidthArgString);
        ss << endl;
        return ss.str();
      }
  };
  
  // ------------------
  //  Parameter packer
  // ------------------
  template<typename T>
    class Param : public ParamName
    {
      private:
        T &value_;
        T const default_value_;

      public:
        using value_type = T;
        Param(std::string description, T& val, 
            std::string short_name, std::string long_name) :
          ParamName(std::move(description), std::move(short_name), std::move(long_name)),
          value_(val), default_value_(val) 
        {}

        T& value() {return value_;}
        T const & value() const {return value_;}
        T const & default_value() const {return default_value_;}
    };

  template<typename T>
  static Param<T> param(char const *description, T& value, char const *short_name, char const *long_name)
  {
    return Param<T>(description,value,short_name,long_name);
  }
  template<typename... Ts>
    Parser<Ts...>
    static pack(int argc, char **argv, Ts&&... ts)
    {
      return Parser<Ts...>(argc, argv, std::forward<Ts>(ts)...);
    }
}
