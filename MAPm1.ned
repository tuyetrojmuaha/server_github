package mm1;

simple MAPs {
	parameters:
	    // XML Description of Arrival Process
		xml model;
		// function to transform the generated interarrival times
		string transform = default("");
		@display("i=block/source");
	gates:
		output out;
}

network MAPm1 {
	submodules:
		maps: MAPs;
		queue: Queue;
		sink: Sink;
	connections:
		maps.out --> queue.in;
		queue.out --> sink.in;
}
