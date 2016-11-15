
package mm1;

simple Queue1 {
	parameters:
		volatile double serviceTime;  // sec
		@display("i=block/queue");
	gates:
		input in;
		output out;
}

simple Queue2 {
	parameters:
		volatile double serviceTime;  // sec
		@display("i=block/queue");
	gates:
		input in;
		output out;
}

simple Queue3 {
	parameters:
		volatile double serviceTime;  // sec
		@display("i=block/queue");
	gates:
		input in;
		output out;
}

simple Queue4 {
	parameters:
		volatile double serviceTime;  // sec
		@display("i=block/queue");
	gates:
		input in;
		output out;
}

network Linear_mm1 {
	submodules:
		gen: Generator;
		queue: Queue;
		queue1: Queue1;
		queue2: Queue2;
		queue3: Queue3;
		queue4: Queue4;
		sink: Sink;
	connections:
		gen.out --> queue.in;
		queue.out --> queue1.in;
		queue1.out --> queue2.in;
		queue2.out --> queue3.in;
		queue3.out --> queue4.in;
		queue4.out --> sink.in;
}