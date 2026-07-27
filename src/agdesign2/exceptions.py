class AgDesignError(Exception):
    """Base exception for antigen design workflow errors."""


class ResolutionError(AgDesignError):
    """Raised when a target cannot be resolved cleanly."""


class ExternalServiceError(AgDesignError):
    """Raised when an external data provider fails."""


class BlastDatabaseError(AgDesignError):
    """Raised when BLAST database setup or execution fails."""
