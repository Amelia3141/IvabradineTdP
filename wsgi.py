import os
import sys
import logging
from api_server import app

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

def initialize_app():
    """Initialize the application and verify the environment"""
    try:
        # Log environment information
        logger.info("Initializing application...")
        logger.info(f"Python version: {sys.version}")
        logger.info(f"Current working directory: {os.getcwd()}")
        
        # Check for required environment variables
        required_vars = ['NCBI_EMAIL', 'NCBI_API_KEY']
        missing_vars = [var for var in required_vars if not os.environ.get(var)]
        
        if missing_vars:
            logger.warning(f"Missing environment variables: {', '.join(missing_vars)}")
        else:
            logger.info("All required environment variables are set")
        
        return app
        
    except Exception as e:
        logger.error(f"Failed to initialize application: {str(e)}")
        raise

if __name__ == "__main__":
    if "--test-only" in sys.argv:
        try:
            initialize_app()
            logger.info("Application initialized successfully in test mode")
            sys.exit(0)
        except Exception as e:
            logger.error(f"Test initialization failed: {str(e)}")
            sys.exit(1)
    else:
        app = initialize_app()
        port = int(os.environ.get('PORT', 10000))
        app.run(host='0.0.0.0', port=port)
