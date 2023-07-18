
#include "../utils/CommandLineArgumentParser.hpp"

#include <iostream>

using namespace std;

int main() {

    cout << endl;
    cout << "Testing the Command Line Argument Parser!!!" << endl;
    cout << endl;

    cout << "Creating the parser object..." << endl;
    CommandLineArgumentParser parser;
    cout << "Created the parser object!" << endl << endl;

    bool successfulOperation;

    cout << "Now, we test the addOption method..." << endl;
    parser.addOption("", "Just a test option 1.", &successfulOperation);
    assert(!successfulOperation);
    parser.addOption("test1", "Just a test option 1.", &successfulOperation);
    assert(successfulOperation);
    parser.addOption("test1", "Just a test option 1.", &successfulOperation);
    assert(!successfulOperation);
    parser.addOption("t2", "Just a test option 2.", &successfulOperation);
    assert(successfulOperation);
    cout << "We have succesfully added two options!" << endl << endl;

    cout << "We print the description of the options..." << endl;
    parser.printOptionDescriptions();
    cout << endl << endl;

    cout << "Now, we parse the commandline arguments..." << endl;
    char *argv1[2] = {"test", "--dne"};
    parser.parseCommandLine(2,  argv1, &successfulOperation);
    assert(!successfulOperation);
    char *argv2[3] = {"test", "--test1", "--test1"};
    parser.parseCommandLine(3,  argv2, &successfulOperation);
    assert(!successfulOperation);
    char *argv3[2] = {"test", "dne"};
    parser.parseCommandLine(2,  argv3, &successfulOperation);
    assert(!successfulOperation);
    parser.resetParser();
    char *argv4[8] = {"test", "--test1", "1", "2", "--help", "--t2", "0.01", "0.5"};
    parser.parseCommandLine(8,  argv4, &successfulOperation);
    assert(successfulOperation);
    cout << endl;

    return 0;

}