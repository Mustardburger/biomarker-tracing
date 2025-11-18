import subprocess, time, re, os, logging
import yaml, copy, argparse

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


def create_yml_files(base_yml: str, diseases: list, save_path: str):
    """
    Create yml files for each disease, then save it to a directory
    """
    os.makedirs(save_path, exist_ok=True)
    with open(base_yml, "r") as f:
        config = yaml.safe_load(f)

    for disease in sorted(diseases):
        new_config = copy.deepcopy(config)
        new_config["inputs"]["disease_name"] = [disease]

        with open(f"{save_path}/{disease}.yml", 'w') as file:
            yaml.dump(new_config, file, default_flow_style=False)


def main(args):
    """
    Run this script if you want to populate yml files across diseases
    """
    # Run code to create yml files
    diseases = os.listdir(args.disease_path)
    diseases = [i.split(".")[0] for i in diseases]
    create_yml_files(args.yml_full_path, diseases, f"{args.save_path}/{args.save_name}")


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--save_path", type=str, help="Save path")
    parser.add_argument("--save_name", type=str, help="Save name")
    parser.add_argument("--yml_full_path", type=str, help="Full path to yml file")
    parser.add_argument("--disease_path", type=str, help="Disease path")
    args = parser.parse_args()
    main(args)
    