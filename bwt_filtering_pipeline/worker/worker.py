#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import sys
import subprocess
import redis
import uuid
import hashlib
import json
import time


class Initializer(object):
    def __init__(self):
        self.queue_name = self.parse_namespace()

    @staticmethod
    def _parse_args():
        p = argparse.ArgumentParser(description="The script exports items from a JSON-based queue to per-node sample data and triggers its processing")
        p.add_argument("-q", "--queue", required=True,
                       help="Redis queue name")
        return p.parse_args()

    def parse_namespace(self):
        n = self._parse_args()
        return n.queue


class Utilities(object):
    @staticmethod
    def ends_with_slash(string):
        return (string + "/", string)[string.endswith("/")]
    @staticmethod
    def get_script_dir():
        return Utilities.ends_with_slash(os.path.dirname(os.path.realpath(sys.argv[0])))
    @staticmethod
    def get_time():
        from datetime import datetime
        now = datetime.now()
        output_list = []
        for time_unit in [now.year, now.month, now.day, now.hour, now.minute, now.second]:
            time_unit = str(time_unit)
            if len(time_unit) < 2:
                time_unit = '0' + time_unit
            output_list.append(time_unit)
        return '-'.join(output_list)
    @staticmethod
    def external_route(args):
        process = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (output, error) = process.communicate()
        process.wait()
        if process.returncode == 0:
            pass
        elif process.returncode == 1:
            print("The command has returned code 1: ", str(args))
        else:
            assert process.returncode > 1
            print("The command has returned code more than 1: ", str(args))
        if error:
            print(error)
        return output.decode("utf-8")


class RedisWQ(object):
    """
    Based on http://peter-hoffmann.com/2012/python-simple-queue-redis-queue.html
    and the suggestion in the redis documentation for RPOPLPUSH, at
    http://redis.io/commands/rpoplpush, which suggests how to implement a work-queue.


    Simple Finite Work Queue with Redis Backend

    This work queue is finite: as long as no more work is added
    after workers start, the workers can detect when the queue
    is completely empty.

    The items in the work queue are assumed to have unique values.

    This object is not intended to be used by multiple threads
    concurrently.
    """
    def __init__(self, name, **redis_kwargs):
        """The default connection parameters are: host='localhost', port=6379, db=0

        The work queue is identified by "name".  The library may create other
        keys with "name" as a prefix.
        """
        self._db = redis.StrictRedis(**redis_kwargs)
        # The session ID will uniquely identify this "worker".
        self._session = str(uuid.uuid4())
        # Work queue is implemented as two queues: main, and processing.
        # Work is initially in main, and moved to processing when a client picks it up.
        self._main_q_key = name
        self._processing_q_key = name + ":processing"
        self._lease_key_prefix = name + ":leased_by_session:"

    def session_id(self):
        """Return the ID for this session."""
        return self._session

    def _main_qsize(self):
        """Return the size of the main queue."""
        return self._db.llen(self._main_q_key)

    def _processing_qsize(self):
        """Return the size of the main queue."""
        return self._db.llen(self._processing_q_key)

    def empty(self):
        """Return True if the queue is empty, including work being done, False otherwise.

        False does not necessarily mean that there is work available to work on right now,
        """
        return self._main_qsize() == 0 and self._processing_qsize() == 0

# TODO: implement this
#    def check_expired_leases(self):
#        """Return to the work queueReturn True if the queue is empty, False otherwise."""
#        # Processing list should not be _too_ long since it is approximately as long
#        # as the number of active and recently active workers.
#        processing = self._db.lrange(self._processing_q_key, 0, -1)
#        for item in processing:
#          # If the lease key is not present for an item (it expired or was
#          # never created because the client crashed before creating it)
#          # then move the item back to the main queue so others can work on it.
#          if not self._lease_exists(item):
#            TODO: transactionally move the key from processing queue to
#            to main queue, while detecting if a new lease is created
#            or if either queue is modified.

    @staticmethod
    def _itemkey(item):
        """Returns a string that uniquely identifies an item (bytes)."""
        return hashlib.sha224(item).hexdigest()

    def _lease_exists(self, item):
        """True if a lease on 'item' exists."""
        return self._db.exists(self._lease_key_prefix + self._itemkey(item))

    def lease(self, lease_secs=60, block=True, timeout=None):
        """Begin working on an item the work queue.
        Lease the item for lease_secs.  After that time, other
        workers may consider this client to have crashed or stalled
        and pick up the item instead.
        If optional args block is true and timeout is None (the default), block
        if necessary until an item is available."""
        if block:
            item = self._db.brpoplpush(self._main_q_key, self._processing_q_key, timeout=timeout)
        else:
            item = self._db.rpoplpush(self._main_q_key, self._processing_q_key)
        if item:
            # Record that we (this session id) are working on a key.  Expire that
            # note after the lease timeout.
            # Note: if we crash at this line of the program, then GC will see no lease
            # for this item a later return it to the main queue.
            itemkey = self._itemkey(item)
            self._db.setex(self._lease_key_prefix + itemkey, lease_secs, self._session)
        return item

    def complete(self, value):
        """Complete working on the item with 'value'.
        If the lease expired, the item may not have completed, and some
        other worker may have picked it up.  There is no indication
        of what happened."""
        self._db.lrem(self._processing_q_key, 0, value)
        # If we crash here, then the GC code will try to move the value, but it will
        # not be here, which is fine.  So this does not need to be a transaction.
        itemkey = self._itemkey(value)
        self._db.delete(self._lease_key_prefix + itemkey, self._session)

# TODO: add functions to clean up all keys associated with "name" when
# processing is complete.
# TODO: add a function to add an item to the queue.  Atomically
# check if the queue is empty and if so fail to add the item
# since other workers might think work is done and be in the process
# of exiting.
# TODO(etune): move to my own github for hosting, e.g. github.com/erictune/rediswq-py and
# make it so it can be pip installed by anyone (see
# http://stackoverflow.com/questions/8247605/configuring-so-that-pip-install-can-work-from-github)
# TODO(etune): finish code to GC expired leases, and call periodically
#  e.g. each time lease times out.


class PipelineWrapper(object):
    def __init__(self, sampledata, refdata, mask, threads, output, no_coverage=False):
        script_dir = Utilities.get_script_dir()
        self._cmd_list = ["python3", '{}nBee.py'.format(script_dir), "-i", sampledata, "-r", refdata, "-m", mask, "-t", threads, "-o", output]
        if no_coverage:
            self._cmd_list = self._cmd_list[:-2] + ["-n"] + self._cmd_list[-2:]

    def run(self):
        print("Started single alignment with arguments: '{}'".format(self._cmd_list))
        Utilities.external_route(self._cmd_list)
        print("Completed single alignment with arguments: '{}'".format(self._cmd_list))


class QueueHandler(object):
    def __init__(self):
        self._queue_name = mainInitializer.queue_name
        self._host_name = "redis"
        self._local_name = subprocess.getoutput("hostname").strip()
        self._queue_connection = RedisWQ(name=self._queue_name, host=self._host_name)
        print("Started worker with sessionID: '{}'".format(self._queue_connection.session_id()))
        print("Initial queue '{a}' state: empty={b}".format(a=self._queue_name, b=str(self._queue_connection.empty())))

    @staticmethod
    def _parse_threads_number(s):
        t = int(subprocess.getoutput("nproc").strip())
        o = 0
        if len(s) > 0:
            if s.isnumeric():
                s = int(s)
                if s > 0:
                    o = min(t, s)
                else:
                    o = t
            else:
                if s == "max":
                    o = t
                elif s == "half":
                    o = int(t / 2)
                elif s == "third":
                    o = int(t / 3)
                elif s == "two_thirds":
                    o = int(t * 2 / 3)
        if o == 0:
            print("Cannot parse the threads number: '{s}'. Using threads number: '{t}'".format(s=s, t=t))
            o = t
        return o

    @staticmethod
    def _dump_sampledata(input_list, output_file):
        output_buffer = ""
        for i in input_list:
            output_buffer += "{a}\t{b}\n".format(a=i["sampledata"]["sample_name"], b=i["sampledata"]["sample_path"])
        with open(file=output_file, mode='r') as f:
            f.write(output_buffer)

    def _run_pipeline(self, json_items_list):
        j = json_items_list[0]
        sampledata_dir = "{a}sampledata_{b}/".format(a=Utilities.ends_with_slash(j["output"]), b=j["mask"])
        os.makedirs(sampledata_dir, exist_ok=True)
        sampledata_file = "{a}{b}_{c}.sampledata".format(a=sampledata_dir, b=self._local_name, c=Utilities.get_time())
        self._dump_sampledata(input_list=json_items_list, output_file=sampledata_file)
        pipeline = PipelineWrapper(sampledata=sampledata_file,
                                   refdata=j[["refdata"]],
                                   mask=j["mask"],
                                   threads=j["threads"],
                                   output=Utilities.ends_with_slash(j["output"]),
                                   no_coverage="no_coverage" in list(j))
        pipeline.run()

    def consume(self):
        json_items_list = []
        threads_number = None
        idle_counter = 0
        max_idle_counter = 100
        while self._queue_connection.empty() and idle_counter < max_idle_counter:
            print("The queue is empty! Paused for 60 seconds, {} attempts left".format(max_idle_counter - idle_counter))
            time.sleep(60)
            idle_counter += 1
            self._queue_connection = RedisWQ(name=self._queue_name, host="redis")
        idle_counter = 0
        while not self._queue_connection.empty() and idle_counter < max_idle_counter:
            queue_item = self._queue_connection.lease(lease_secs=10, block=True, timeout=2)
            if queue_item is not None:
                string_item = queue_item.decode("utf=8")
                try:
                    json_item = json.loads(string_item)
                    # Note that all entries must have the same refdata
                    # JSON keys: {"refdata": "", "sampledata": {"sample_name": "", "sample_path": ""}, "mask": "", "threads": "", "output": ""}
                    json_items_list.append(json_item)
                    if len(json_items_list) == 1:
                        threads_number = self._parse_threads_number(json_item["threads"])
                    # A pause to allow other nodes access the queue
                    time.sleep(5)
                    if threads_number and len(json_items_list) == threads_number:
                        print("Loaded full queue on '{}'".format(self._local_name))
                        self._run_pipeline(json_items_list)
                        print("Processing is finished, loading next queue on '{}'".format(self._local_name))
                        json_items_list = []
                except ValueError:
                    print("Cannot parse JSON: '{}'".format(string_item))
                self._queue_connection.complete(queue_item)
                idle_counter = 0
            else:
                print("Waiting for work, {} attempts left".format(max_idle_counter - idle_counter))
                idle_counter += 1
        if len(json_items_list) > 0:
            print("Processing the last queue on '{}'".format(self._local_name))
            self._run_pipeline(json_items_list)
            print("Queue is empty, exiting")


if __name__ == '__main__':
    # The script part is based on: https://raw.githubusercontent.com/kubernetes/website/master/docs/tasks/job/fine-parallel-processing-work-queue/worker.py
    # Uncomment next two lines if you do not have Kube-DNS working.
    # import os
    # host = os.getenv("REDIS_SERVICE_HOST")

    mainInitializer = Initializer()
    handler = QueueHandler()
    handler.consume()
