#include <iostream>
#include <list>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <algorithm> // for std::copy

#import "C:\Program Files (x86)\Agilent Technologies\ChemStation\CORE\AgtExternalIfcBroker.exe" rename("GetObject","_GetObject") no_namespace named_guids raw_interfaces_only

#define CHECK(hresult, ErrOut, retCode) if (FAILED(hresult)){std::cout << ErrOut; return retCode;}

int Connect(int PCSNum, bool online, IPublishedChemStation** pChemStation)
{
    *pChemStation = NULL;

    // get the external interface broker
    IPublishedChemStationsPtr pPublishedChemStations;
    CHECK(pPublishedChemStations.CreateInstance(__uuidof(ChemStationBroker), NULL), "Fatal: Could not instantiate broker.\n", 43);

    std::cout << "Got the external interface broker.\n";
    long count;
    CHECK(pPublishedChemStations->get_Count(&count), "Fatal: Automation Collection did not return count.\n", 44);
    std::cout << "Found " << count << " running ChemStation instance(s).\n";

    //looking for requested ChemStation instance
    VARIANT idx;
    idx.vt = VT_I4;
    for (idx.intVal = 1; idx.intVal <= count; idx.intVal++)
    {
        IDispatchPtr Item;
        CHECK((pPublishedChemStations->Item(idx, &Item)), "Fatal: Automation Collection did not return item for valid index.\n", 45);

        IChemStationConnectorPtr csc;
        CHECK(Item.QueryInterface(__uuidof(IChemStationConnector), &csc), "Fatal: Automation Collection item did not support expected interface.\n", 46);

        long i;
        CHECK(csc->get_PCSNum(&i), "Fatal: Chemstation Connector did not expose PCS number.\n", 47);

        VARIANT_BOOL o;
        CHECK(csc->get_Online(&o), "Fatal: Chemstation Connector did not expose online flag.\n", 48);

        if ((i == PCSNum) && (online == (o == VARIANT_TRUE)))
        {
            std::cout << "Found requested ChemStation instance, authenticating\n";
            // Authenticate
            CHECK(csc->Connect2(_bstr_t(""), _bstr_t(""), _bstr_t(""), pChemStation), "Failed to connect.\n", 49);
            std::cout << "Connected to requested ChemStation instance.\n";

            return 0;
        }
    }
    std::cout << "Did not find requested ChemStation instance.\n";
    return 50;
}

int main(int argc, char** argv)
{

    // Read txt and make VECTOR

    std::ifstream is("D:\\Jim_Boelrijk\\AutoLC\\AutoLC-BO\\parameters_4step.txt");



    std::istream_iterator<double> start(is), end;
    std::vector<double> numbers(start, end);
    std::cout << "Read " << numbers.size() << " numbers" << std::endl;

    // print the numbers to stdout
   // std::cout << "numbers read in:\n";
   // std::copy(numbers.begin(), numbers.end(),
   //     std::ostream_iterator<double>(std::cout, " "));
   // std::cout << std::endl;
  //  std::cout << numbers[0]  << std::endl;



// Make connection
    if (FAILED(CoInitializeEx(NULL, COINIT_APARTMENTTHREADED)))
    {
        std::cout << "Fatal: Could not initialize COM.\n";
        return 42;
    }

    int PCSNum = 1;
    bool online = true;
    if (argc > 1)
        sscanf_s(argv[1], "%d", &PCSNum);

    if (argc > 2)
    {
        int i;
        sscanf_s(argv[2], "%d", &i);
        online = (i != 0);
    }
    std::cout << "Trying to connect to ChemStation " << PCSNum << ((online) ? " online\n" : " offline\n");

    // excecute commands
    IPublishedChemStationPtr CS;
    int retCode = Connect(PCSNum, online, &CS);
    if (retCode == 0)
    {
        IChemStationCPPtr CP;
        CHECK(CS->get_CP(&CP), "Failed to retrieve access to command processor.\n", 51);
        std::cout << "Got access to command processor.\n";

        //TODO: Implement all these lines for the desired settings
        //TODO: This needs be done differently for the first line of timetable, we need to have a look what the variable names are for this.
        //TODO: main() needs to accept a list of LC settings which we can then place as variables into the strings below

        //load macros
        //CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t("macro \"EXPORTSIGNAL_ELSD.mac\"")), "Failed to load hook macro.\n", 52);
        //CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t("macro \"Export3D_seq.mac\"")), "Failed to load hook macro.\n", 52); //Define path within the macro!
        CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t("macro \"Export3D_BO.mac\"")), "Failed to load hook macro.\n", 52); //Define path within the macro!

        CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t("macro \"WhosTurn_BO.mac\"")), "Failed to load hook macro.\n", 52);
        //CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t("SetHook \"PostRun\", \"EXPORTSIGNAL_ELSD\"")), "Failed to load post run hook.\n", 52);

        // TODO: Make sure blank and real measurement are written to unique file.
        // 
        //CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t("SetHook \"PostRun\", \"Export3D\", \"D:\\Jim_Boelrijk\\AutoLC\\AutoLC-BO\\Export3D.csv\"")), "Failed to load export 3d hook.\n", 52);

        CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t("SetHook \"PostRun\", \"Export3D_BO\"")), "Failed to load export 3d hook.\n", 52);
        //CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t("SetHook \"PostRun\", \"Export3D_Jim\"")), "Failed to load export 3d hook.\n", 52);

        CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t("SetHook \"PostSeq\", \"WhosTurn_BO\"")), "Failed to load export waiter hook.\n", 52);

        CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t("AUTOMETHOD = 1")), "Failed to load hook.\n", 52);
        // Make sequence in commands

        //Export3D "O:\Boelrijk_Jim\AutoLC\Export3D.csv"

        // TODO: change method names
        //CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t("LoadSequence \"C:\\Users\\Public\\Documents\\ChemStation\\1\\Sequence\\\", \"ELSD_AUTOLC.S\"")), "Failed load ELSD_AUTOLC.S sequence.\n", 57);
        //CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t("LoadMethod \"C:\\Users\\Public\\Documents\\ChemStation\\1\\Methods\\\", \"AutoELSD.M\"")), "Failed loading AutoELSD.M method.\n", 57);

        CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t("LoadSequence \"C:\\Users\\Public\\Documents\\ChemStation\\1\\Sequence\\\", \"BO-Seq.S\"")), "Failed load Opt sequence.\n", 57);
        CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t("LoadMethod \"C:\\Users\\Public\\Documents\\ChemStation\\1\\Methods\\\", \"Auto-LC-BO.M\"")), "Failed load method.\n", 57);

        // set variables
        //CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t(("SetObjHdrVal RCPMP1Method[1],\"Flow\"," + std::to_string(numbers[0])).c_str())), "Failed to load flow.\n", 53); // set flow through input
        //CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t(("SetObjHdrVal RCPMP1Method[1],\"Flow\"," + std::to_string(0.65)).c_str())), "Failed to load flow.\n", 53); //set flow manually

       //CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t(("SetObjHdrVal RCPMP1Method[1],\"PumpChannel2_CompositionPercentage\"," + std::to_string(numbers[0])).c_str())), "Failed to load phi init.\n", 54);
        //CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t(("SetTabVal  RCPMP1Method[1],\"TimeTable\", 1, \"SolventCompositionPumpChannel2_Percentage\"," + std::to_string(numbers[0])).c_str())), "Failed to load phi init2.\n", 55);

        // Moved up all indices, because index 1 now is only used to activate the flow rate, index 0 sets idle flow.
        CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t(("SetObjHdrVal RCPMP1Method[1],\"PumpChannel2_CompositionPercentage\"," + std::to_string(0)).c_str())), "Failed to load phi init.\n", 54); //TODO: fix phi init
        CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t(("SetTabVal  RCPMP1Method[1],\"TimeTable\", 2, \"SolventCompositionPumpChannel2_Percentage\"," + std::to_string(0)).c_str())), "Failed to load phi init2.\n", 55); //TODO: fix phi init
        CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t(("SetTabVal  RCPMP1Method[1],\"TimeTable\", 2, \"Time\"," + std::to_string(0.25)).c_str())), "Failed to load init time.\n", 56);
        //CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t(("SetTabVal  RCPMP1Method[1],\"TimeTable\", 2, \"FlowFlow\"," + std::to_string(0.65)).c_str())), "Failed to load flow 1.\n", 56);

        CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t(("SetTabVal  RCPMP1Method[1],\"TimeTable\", 3, \"SolventCompositionPumpChannel2_Percentage\"," + std::to_string(numbers[0])).c_str())), "Failed to load phi1 2.\n", 55);
        CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t(("SetTabVal  RCPMP1Method[1],\"TimeTable\", 3, \"Time\"," + std::to_string(numbers[3] + 0.25)).c_str())), "Failed to load time step 1.\n", 56);
        //CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t(("SetTabVal  RCPMP1Method[1],\"TimeTable\", 3, \"FlowFlow\"," + std::to_string(0.65)).c_str())), "Failed to load flow 2.\n", 56);

        //CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t(("SetTabVal  RCPMP1Method[1],\"TimeTable\", 3, \"SolventCompositionPumpChannel2_Percentage\"," + std::to_string(numbers[3])).c_str())), "Failed to load phi final 2.\n", 55);
        CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t(("SetTabVal  RCPMP1Method[1],\"TimeTable\", 4, \"SolventCompositionPumpChannel2_Percentage\"," + std::to_string(numbers[1])).c_str())), "Failed to load phi2.\n", 55); //TODO: fix phi final
        CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t(("SetTabVal  RCPMP1Method[1],\"TimeTable\", 4, \"Time\"," + std::to_string(0.25 + numbers[4])).c_str())), "Failed to load time step 2\n", 56);
        //CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t(("SetTabVal  RCPMP1Method[1],\"TimeTable\", 4, \"FlowFlow\"," + std::to_string(0.65)).c_str())), "Failed to load flow 3.\n", 56);

        CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t(("SetTabVal  RCPMP1Method[1],\"TimeTable\", 5, \"SolventCompositionPumpChannel2_Percentage\"," + std::to_string(numbers[2])).c_str())), "Failed to load phi3.\n", 55); //TODO: fix phi final
        CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t(("SetTabVal  RCPMP1Method[1],\"TimeTable\", 5, \"Time\"," + std::to_string(0.25 + numbers[5])).c_str())), "Failed to load time step 3\n", 56);
        //CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t(("SetTabVal  RCPMP1Method[1],\"TimeTable\", 5, \"FlowFlow\"," + std::to_string(0.65)).c_str())), "Failed to load flow 4.\n", 56);

        CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t(("SetTabVal  RCPMP1Method[1],\"TimeTable\", 6, \"SolventCompositionPumpChannel2_Percentage\"," + std::to_string(100)).c_str())), "Failed to load phi recalibration 2.\n", 55);
        CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t(("SetTabVal  RCPMP1Method[1],\"TimeTable\", 6, \"Time\"," + std::to_string(20 + 0.25)).c_str())), "Failed to load recalibration time.\n", 56);
        //CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t(("SetTabVal  RCPMP1Method[1],\"TimeTable\", 6, \"FlowFlow\"," + std::to_string(0.65)).c_str())), "Failed to load flow 5.\n", 56);


        CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t(("SetTabVal  RCPMP1Method[1],\"TimeTable\", 7, \"SolventCompositionPumpChannel2_Percentage\"," + std::to_string(100)).c_str())), "Failed to load phi recalibration 3.\n", 55);
        CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t(("SetTabVal  RCPMP1Method[1],\"TimeTable\", 7, \"Time\"," + std::to_string(20 + 0.25 + 0.1)).c_str())), "Failed to load recalibration time.\n", 56);

        CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t(("SetTabVal  RCPMP1Method[1],\"TimeTable\", 8, \"SolventCompositionPumpChannel2_Percentage\"," + std::to_string(0)).c_str())), "Failed to load phi recalibration 4.\n", 55);
        CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t(("SetTabVal  RCPMP1Method[1],\"TimeTable\", 8, \"Time\"," + std::to_string(20 + 0.25 + 0.1 + 1.5)).c_str())), "Failed to load recalibration time.\n", 56);
        //CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t(("SetTabVal  RCPMP1Method[1],\"TimeTable\", 7, \"FlowFlow\"," + std::to_string(0.01)).c_str())), "Failed to load flow 6.\n", 56);

        // TODO: add another line where we go back to 0 mobile phase.

        // Set the stop time of the experiment to the last gradient time + 2 extra minutes
        CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t(("SetObjHdrVal RCPMP1Method[1],\"StopTime_Time\"," + std::to_string(20 + 0.25 + 3)).c_str())), "Failed to set stop time.\n", 54);

        // upload and start method
        CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t("DownloadRCMethod PMP1")), "Failed download to LC.\n", 57);

        CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t("SaveMethod \"C:\\Users\\Public\\Documents\\ChemStation\\1\\Methods\\\", \"Auto-LC-BO.M\", \" \" ")), "Failed save method.\n", 57);
        CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t("StartSequence")), "Failed start measurent.\n", 57);
        //CHECK(CP->RequestWait(_bstr_t("CP"), _bstr_t("StartMethod")), "Failed to start.\n", 58);

        std::cout << "Success: Executed CP command string.\n";

    }

    return retCode;
}

