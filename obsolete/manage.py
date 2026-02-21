import subprocess
import sys
import os
import signal
import time
import argparse
import socket

def is_port_in_use(port):
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        return s.connect_ex(('localhost', port)) == 0

def get_status():
    print("\n--- Service Status ---")
    backend_running = is_port_in_use(8000)
    frontend_running = is_port_in_use(3000)
    
    print(f"Backend (Port 8000): {'[RUNNING]' if backend_running else '[STOPPED]'}")
    print(f"Frontend (Port 3000): {'[RUNNING]' if frontend_running else '[STOPPED]'}")
    print("----------------------\n")
    return backend_running, frontend_running

def start_backend():
    if is_port_in_use(8000):
        print("Backend is already running on port 8000.")
        return None
    print("Starting backend...")
    # Using uv to run if available, otherwise fallback to sys.executable
    cmd = ["uv", "run", "python", "-m", "backend.main"]
    try:
        proc = subprocess.Popen(
            cmd,
            stdout=open("backend.log", "a"),
            stderr=open("backend.log", "a"),
            preexec_fn=os.setsid
        )
        print("Backend process initiated. PID:", proc.pid)
        return proc
    except FileNotFoundError:
        print("uv not found, falling back to standard python...")
        proc = subprocess.Popen(
            [sys.executable, "-m", "backend.main"],
            stdout=open("backend.log", "a"),
            stderr=open("backend.log", "a"),
            preexec_fn=os.setsid
        )
        return proc

def start_frontend():
    if is_port_in_use(3000):
        print("Frontend is already running on port 3000.")
        return None
    print("Starting frontend (dev mode)...")
    proc = subprocess.Popen(
        ["npm", "start"],
        cwd="frontend",
        env={**os.environ, "BROWSER": "none"},
        stdout=open("frontend.log", "a"),
        stderr=open("frontend.log", "a"),
        preexec_fn=os.setsid
    )
    print("Frontend process initiated. PID:", proc.pid)
    return proc

def stop_all():
    print("Stopping all services...")
    # Kill local processes
    subprocess.run(["pkill", "-f", "backend.main"], stderr=subprocess.DEVNULL)
    subprocess.run(["pkill", "-f", "react-scripts"], stderr=subprocess.DEVNULL)
    
    # Kill docker if running
    if os.path.exists("docker-compose.yml"):
        print("Stopping Docker containers...")
        subprocess.run(["docker-compose", "down"], stderr=subprocess.DEVNULL)
    
    print("Stop signal sent to backend, frontend, and Docker processes.")

def rebuild_data():
    print("\n--- Rebuilding and Cleaning Data ---")
    
    # 1. Run process_data.py
    print("Step 1: Running process_data.py...")
    subprocess.run(["uv", "run", "--with", "pandas", "python", "scripts/process_data.py"], check=True)
    
    # 2. Run finalize_master.py (this includes our new cleaning logic)
    print("Step 2: Running finalize_master.py (with gene/OMIM cleaning)...")
    subprocess.run(["uv", "run", "--with", "pandas", "python", "scripts/finalize_master.py"], check=True)
    
    # 3. Run generate_phenopackets.py
    print("Step 3: Generating Phenopackets JSON files...")
    subprocess.run(["uv", "run", "--with", "pandas", "python", "scripts/generate_phenopackets.py"], check=True)
    
    # 4. Reset database
    if os.path.exists("phenopackets.db"):
        print("Step 4: Resetting database (deleting phenopackets.db)...")
        os.remove("phenopackets.db")
    else:
        print("Step 4: Database file not found, skipping deletion.")
        
    print("Data rebuild complete.\n")

def main():
    parser = argparse.ArgumentParser(description="Manage Saudi Phenopackets Search App")
    parser.add_argument("command", choices=["start", "stop", "restart", "logs", "status", "rebuild", "docker-start"], help="Command to execute")
    
    args = parser.parse_args()

    if args.command == "start":
        start_backend()
        start_frontend()
        print("\nStartup sequence complete. It may take 15-20s for the HPO ontology to load.")
        print("Check status with: python manage.py status")
    elif args.command == "stop":
        stop_all()
    elif args.command == "restart":
        stop_all()
        time.sleep(2)
        start_backend()
        start_frontend()
        print("Services restarted.")
    elif args.command == "status":
        get_status()
    elif args.command == "logs":
        print("Following logs (Ctrl+C to stop)...")
        try:
            subprocess.run(["tail", "-f", "backend.log", "frontend.log"])
        except KeyboardInterrupt:
            pass
    elif args.command == "rebuild":
        stop_all()
        rebuild_data()
        start_backend()
        start_frontend()
        print("Data rebuilt and services restarted.")
    elif args.command == "docker-start":
        if os.path.exists("docker-compose.yml"):
            print("Starting services with Docker Compose...")
            subprocess.run(["docker-compose", "up", "--build", "-d"])
        else:
            print("docker-compose.yml not found.")

if __name__ == "__main__":
    main()
