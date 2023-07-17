
#include "CommandLineArgumentParser.hpp"

#include <iostream>

CommandLineArgumentParser::CommandLineArgumentParser() {
    Option* help = new Option;
    help -> isSet = false;
    help -> num_arguments = -1;
    help -> description = "Prints the preamble and options list for the program.";
    arguments["--help"] = help;
}

string form_option(string option) {
    if (option.length() == 1) {
        return "-" + option;
    }
    return "--" + option;
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

template <class T>
T* CommandLineArgumentParser::getOptionArguments(string option, int* num_arguments, bool* successfulOperation) {

    option = form_option(option);

    if (!arguments.count(option)) {
        cerr << "Option " << option << " does not exists!" << endl;
        *successfulOperation = false;
        return;
    }

    Option* opt = arguments[option];

    if (!opt -> isSet || opt -> num_arguments == 0) {
        *successfulOperation = true;
        return NULL;
    }

    T* args = new T[opt -> num_arguments];

    string temp;

    for (int i = 0; i < opt -> num_arguments; i++) {
        try {
            temp = opt -> option_arguments.front();
            opt -> option_arguments.pop_front();
            temp >> args[i];
            opt -> option_arguments.push_back(temp);
        } catch (...) {
            cerr << "Invalid type cast when getting arguments for " << option << " option." << endl;
            *successfulOperation = false;
        }
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

void CommandLineArgumentParser::resetParser() {

    for (map<string, Option*>::iterator it = arguments.begin(); it != arguments.end(); it++) {
        it -> second -> isSet = false;
        for (int i = 0; i < it -> second -> num_arguments; i++) {
            it -> second -> option_arguments.pop_front();
        }
        it -> second -> num_arguments = -1;
    }
    
}

CommandLineArgumentParser::~CommandLineArgumentParser() {

    for (map<string, Option*>::iterator it = arguments.begin(); it != arguments.end(); it++) {
        delete it -> second;
    }

}