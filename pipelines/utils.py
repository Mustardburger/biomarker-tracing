import subprocess, time, re, os, logging

def submit_job_and_wait(bsub_cmd: list, wait_time=10):
    # Submit job and capture job ID
    result = subprocess.run(bsub_cmd, capture_output=True, text=True)
    m = re.search(r"<(\d+)>", result.stdout)
    jobid = m.group(1)
    logging.error(f"Submitted job {jobid}, waiting...")

    # Sleep for about 10s to wait for the job to get registered
    time.sleep(10)

    # Poll until job finishes
    while True:
        bjobs = subprocess.run(["bjobs", jobid], capture_output=True, text=True)
        lines = bjobs.stdout.strip().splitlines()
        fields = lines[1].split()

        # 0=JOBID, 1=USER, 2=JOB_NAME, 3=STAT
        stat = fields[3]
        if stat == "DONE" or stat == "EXIT":
            logging.error(f"stat={stat}, job {jobid} finished")
            break
        else:
            logging.error(f"stat={stat}, job {jobid} still running or pending")
        time.sleep(wait_time)