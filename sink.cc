#include <string.h>
#include <omnetpp.h>
class Sink: public cSimpleModule{
private:
    // online stats
    cStdDev delayStats;
    cOutVector delayVector;

    simtime_t a=0;
    cOutVector deparVector;

    cOutVector line_delay;

public:
    Sink();
    virtual~ Sink();
protected:
    virtual void initialize();
    virtual void finish();
    virtual void handleMessage(cMessage *msg);
};
Define_Module(Sink);
Sink::Sink(){
}
Sink::~Sink(){
}
void Sink::initialize(){
    delayStats.setName("TotalDelay");
    delayVector.setName("Delay");

    deparVector.setName("departure_time");
    line_delay.setName("line_delay");
}
void Sink::finish(){
    recordScalar("Ave_delay",delayStats.getMean());
    recordScalar("Number_of_packet",delayStats.getCount());
}
void Sink::handleMessage(cMessage *msg){
    // compute queueing delay
    simtime_t delay = simTime() - msg -> getCreationTime();
    // update stats
    delayStats.collect(delay);
    delayVector.record(delay);


    deparVector.record(simTime()-a);
    a = simTime();

    line_delay.record(delayStats.getMean());

    //delete msg
    delete(msg);
}
