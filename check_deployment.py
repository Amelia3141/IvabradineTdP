import requests
import time
import sys

def check_deployment(url, max_retries=10, delay=30):
    """Check if the deployment is successful by hitting the health endpoint"""
    print(f"Checking deployment at {url}")
    
    for i in range(max_retries):
        try:
            print(f"\nAttempt {i + 1}/{max_retries}...")
            response = requests.get(f"{url}/health")
            
            if response.status_code == 200:
                data = response.json()
                print("\nDeployment successful! Service health check response:")
                print(f"Status: {data.get('status')}")
                print(f"Components: {data.get('components', {})}")
                print(f"Environment: {data.get('environment', {})}")
                print(f"Versions: {data.get('versions', {})}")
                return True
                
        except requests.RequestException as e:
            print(f"Request failed: {e}")
        
        if i < max_retries - 1:
            print(f"Retrying in {delay} seconds...")
            time.sleep(delay)
    
    print("\nDeployment check failed after maximum retries")
    return False

if __name__ == "__main__":
    url = "https://ivabradinetdp.onrender.com"
    success = check_deployment(url)
    sys.exit(0 if success else 1)
