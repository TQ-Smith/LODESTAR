
#ifndef _COMMAND_LINE_ARGUMENT_PARSER_HPP_
#define _COMMAND_LINE_ARGUMENT_PARSER_HPP_

#include <string>
#include <sstream>
#include <list>
#include <map>
#include <iostream>

using namespace std;

class CommandLineArgumentParser {

    public:

        CommandLineArgumentParser();

        void addOption(string option, string description, bool* successfulOperation);

        template <typename T>
        T* getOptionArguments(string option, int* num_arguments, bool* successfulOperation);

        void parseCommandLine(int argc, char *argv[], bool* successfulOperation);

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

string form_option(string option) {
    if (option.length() == 1) {
        return "-" + option;
    }
    return "--" + option;
}

CommandLineArgumentParser::CommandLineArgumentParser() {
    bool successfulOperation;
    addOption("help",  "Prints the preamble and options list for the program.", &successfulOperation);
}

void CommandLineArgumentParser::addOption(string option, string description, bool* successfulOperation) {

    if (option.length() == 0) {
        cerr << "Cannot use the empty string as an option!" << endl;
        *successfulOperation = false;
        return;
    }

    option = form_option(option);

    if (arguments.count(option)) {
        cerr << "Option " << option << " already exists!" << endl;
        *successfulOperation = false;
        return;
    }

    Option* temp = new Option;
    temp -> isSet = false;
    temp -> num_arguments = -1;
    temp -> description = description;

    arguments[option] = temp;

    *successfulOperation = true;

}

void CommandLineArgumentParser::parseCommandLine(int argc, char *argv[], bool* successfulOperation) {

    Option* opt = NULL;

    string option;

    for (int i = 1; i < argc; i++) {

        if (argv[i][0] == '-') {

            option = argv[i];

            if (!arguments.count(option)) {
                cerr << "Option " << option << " does not exists!" << endl;
                *successfulOperation = false;
                return;
            }

            opt = arguments[option];

            if (opt -> isSet) {
                cerr << "Option " << option << " was already used!" << endl;
                *successfulOperation = false;
                return;
            }

            opt -> isSet = true;
            opt -> num_arguments = 0;

        } else {

            if (opt == NULL) {
                cerr << "No option created for " << argv[i] << "." << endl;
                *successfulOperation = false;
                return;
            }

            opt -> option_arguments.push_back(string(argv[i]));
            opt -> num_arguments++;

        }

    }

    *successfulOperation = true;

}

template <typename T>
T* CommandLineArgumentParser::getOptionArguments(string option, int* num_arguments, bool* successfulOperation) {

    if (!arguments.count(option)) {
        cerr << "Option " << option << " does not exists!" << endl;
        *successfulOperation = false;
        return NULL;
    }

    Option* opt = arguments[option];

    if (!opt -> isSet || opt -> num_arguments == 0) {
        *num_arguments = opt -> num_arguments;
        *successfulOperation = true;
        return NULL;
    }

    T* args = new T[opt -> num_arguments];
    
    int index = 0;
    for (list<string>::iterator it = opt -> option_arguments.begin(); it != opt -> option_arguments.end(); it++) {
        istringstream strstream_arg(*it);
        if (!(strstream_arg >> args[index])) {
            cerr << "Cannot perform type cast of argument for option " << option << "!" << endl;
            delete [] args;
            *successfulOperation = false;
            return NULL;
        }
        index++;
    }

    *num_arguments = opt -> num_arguments;

    *successfulOperation = true;

    return args;

}

void CommandLineArgumentParser::printOptionDescriptions() {

    for (map<string, Option*>::iterator it = arguments.begin(); it != arguments.end(); it++) {
        cout << it -> first << ": " << it -> second -> description << endl;
    }

}

CommandLineArgumentParser::~CommandLineArgumentParser() {

    for (map<string, Option*>::iterator it = arguments.begin(); it != arguments.end(); it++) {
        delete (it -> second);
    }

}

#endif