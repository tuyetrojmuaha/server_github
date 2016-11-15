package mm1;

simple Generator {
	parameters:
		volatile double interArrivalTime; // sec
		@display("i=block/source");
	gates:
		output out;
}

simple Queue {
	parameters:
		volatile double serviceTime;  // sec
		@display("i=block/queue");
	gates:
		input in;
		output out;
}

simple Sink {
	parameters:
		@display("i=block/sink");
	gates:
		input in;
}

network MM1 {
	submodules:
		gen: Generator;
		queue: Queue;
		sink: Sink;
	connections:
		gen.out --> queue.in;
		queue.out --> sink.in;
}