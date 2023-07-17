
#ifndef _COMMAND_LINE_ARGUMENT_PARSER_HPP_
#define _COMMAND_LINE_ARGUMENT_PARSER_HPP_

#include <string>
#include <list>
#include <map>

using namespace std;

class CommandLineArgumentParser {

    public:

        CommandLineArgumentParser();

        void addOption(string option, string description, bool* successfulOperation);

        template <class T>
        T* getOptionArguments(string option, int* num_arguments, bool* successfulOperation);

        void parseCommandLine(int argc, char *argv[], bool* successfulOperation);

        void resetParser();

        void printOptionDescriptions();

        ~CommandLineArgumentParser();

    private:

        struct Option {
            bool isSet;
            int num_arguments;
            string description;
            list<string> option_arguments;
        };

        map<string, Option*> arguments;

};

#endif