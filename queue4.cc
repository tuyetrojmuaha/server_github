#include <string.h>
#include <omnetpp.h>
class Queue4 : public cSimpleModule {
private:
    cQueue buffer;
    cMessage *endServiceEvent;
    simtime_t service_time;

    cStdDev waitStats;
    cOutVector waitVector;

    cStdDev lengthStats;
    cOutVector lengthVector;
public:
    Queue4();
    virtual~ Queue4();
protected:
    virtual void initialize();
    virtual void finish();
    virtual void handleMessage(cMessage *msg);
};
Define_Module(Queue4);
Queue4::Queue4(){
    endServiceEvent = NULL;
}
Queue4::~Queue4(){
    cancelAndDelete(endServiceEvent);
}
void Queue4::initialize(){
    endServiceEvent = new cMessage("endService");

    lengthStats.setName("Total");
    lengthVector.setName("queue_length");

    waitStats.setName("Total_w");
    waitVector.setName("wait_length");
}
void Queue4::finish(){
    recordScalar("Ave_queue_length",lengthStats.getMean());
    recordScalar("Number_of_packet_in_queue",lengthStats.getCount());

    recordScalar("Ave_wait_length",waitStats.getMean());
    recordScalar("Number_of_packet_in_wait",waitStats.getCount());
}
void Queue4::handleMessage(cMessage *msg){
    cMessage *pkt;
    // if msg is endServiceEvent, then
    // dequeue and send the pkt to the output
    // if another pkt is available in the buffer, then
    // start a new service
    // if msg is a packet, then
    // enqueue the pkt
    // if server idling, then
    // start a new service
    if (msg == endServiceEvent){
        // dequeue
        pkt = (cMessage*)buffer.pop();
        //send
        send(pkt,"out");
        if (!buffer.empty()) {
            //if another pkt is available
            // start the service
            service_time = par("serviceTime");
            scheduleAt(simTime()+service_time,endServiceEvent);
        }
    } else {

        lengthStats.collect(buffer.length());
        lengthVector.record(lengthStats.getMean());

        // msg is a packet enqueue
        buffer.insert(msg);
        //if the service is idling
        if (!endServiceEvent ->isScheduled()) {

            waitStats.collect(0);
            waitVector.record(waitStats.getMean());

            //start the service
            service_time = par("serviceTime");
            scheduleAt(simTime()+service_time,endServiceEvent);


        }
        else {
            waitStats.collect(buffer.length()-2);
            waitVector.record(waitStats.getMean());
        }
    }

}
