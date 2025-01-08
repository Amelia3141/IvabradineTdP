import requests
import time
import sys
import json

def check_deployment(url, max_retries=20, initial_delay=60, subsequent_delay=30):
    """Check if the deployment is successful by hitting various endpoints"""
    print(f"\nChecking deployment at {url}")
    print(f"Initial delay: {initial_delay}s, Subsequent delays: {subsequent_delay}s")
    print(f"Maximum retries: {max_retries}")
    
    # Initial delay to allow for startup
    print(f"\nWaiting {initial_delay}s for initial deployment...")
    time.sleep(initial_delay)
    
    endpoints = {
        "health": "/health",
        "root": "/"
    }
    
    for i in range(max_retries):
        print(f"\nAttempt {i + 1}/{max_retries}...")
        success = True
        
        for name, path in endpoints.items():
            try:
                print(f"\nChecking {name} endpoint ({url}{path})...")
                response = requests.get(f"{url}{path}", timeout=30)
                
                print(f"Status code: {response.status_code}")
                if response.status_code == 200:
                    data = response.json()
                    print(f"Response data: {json.dumps(data, indent=2)}")
                    
                    if name == "health":
                        # Verify key health components
                        status = data.get('status')
                        components = data.get('components', {})
                        env = data.get('environment', {})
                        
                        print("\nHealth check details:")
                        print(f"Status: {status}")
                        print(f"Components: {components}")
                        print(f"Environment: {env}")
                        
                        if not all([env.get('NCBI_EMAIL'), env.get('NCBI_API_KEY')]):
                            print("Warning: Missing required environment variables")
                            success = False
                else:
                    print(f"Error: Unexpected status code {response.status_code}")
                    success = False
                    
            except requests.RequestException as e:
                print(f"Request failed: {e}")
                success = False
                continue
        
        if success:
            print("\n✅ Deployment check successful! All endpoints are responding correctly.")
            return True
            
        if i < max_retries - 1:
            print(f"\nSome checks failed. Retrying in {subsequent_delay} seconds...")
            time.sleep(subsequent_delay)
    
    print("\n❌ Deployment check failed after maximum retries")
    return False

if __name__ == "__main__":
    url = "https://ivabradinetdp.onrender.com"
    success = check_deployment(url)
    sys.exit(0 if success else 1)
