package mm1;

network Linear_MAPm1 {
	submodules:
		maps: MAPs;
		queue: Queue;
		queue1: Queue1;
		queue2: Queue2;
		queue3: Queue3;
		queue4: Queue4;
		sink: Sink;
	connections:
		maps.out --> queue.in;
		queue.out --> queue1.in;
		queue1.out --> queue2.in;
		queue2.out --> queue3.in;
		queue3.out --> queue4.in;
		queue4.out --> sink.in;
}