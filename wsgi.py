from api_server import app
import sys

if __name__ == "__main__":
    if "--test-only" in sys.argv:
        print("API server initialized successfully")
        sys.exit(0)
    app.run()
